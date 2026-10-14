"""Alternate NumPy-only equations; no imports from the project implementation.

Inputs: frozen calendar, receiver policies, calibration counts and references.
Reconstructs the full source and propagates the complete eight-qubit contrast.
This is an internal alternate implementation, not independent external review.
"""
import base64
import io
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)

def kron(items):
    out=np.ones((1,1),complex)
    for item in items:out=np.kron(out,item)
    return out

def source(seed):
    gamma=[kron([Z]*q+[p]+[I]*(3-q)) for q in range(4) for p in (X,Y)]
    parity=kron([Z]*4)
    draws=np.random.default_rng(seed).normal(0,math.sqrt(6/8**3),70)
    h=np.zeros((16,16),complex)
    for indices,j in zip(itertools.combinations(range(8),4),draws):
        g=np.eye(16,dtype=complex)
        for q in indices:g=g@gamma[q]
        h-=j*g/4
    values,vectors=np.linalg.eigh(h)
    u=(vectors*np.exp(-1j*values))@vectors.conj().T
    thermal_filter=(vectors*np.exp(-2*(values-values.min())))@vectors.conj().T
    interaction=sum(.5j*np.kron(g@parity,g) for g in gamma)
    _,basis=np.linalg.eigh(interaction)
    vector=(thermal_filter@basis[:,0].reshape(16,16)).reshape(256)
    vector/=np.linalg.norm(vector)
    tfd=np.outer(vector,vector.conj())
    thermal_weights=np.exp(-4*(values-values.min()))
    thermal=(vectors*(thermal_weights/thermal_weights.sum()))@vectors.conj().T
    thermal_errors=[float(np.max(abs(np.trace(tfd.reshape(16,16,16,16),axis1=a,axis2=b)-thermal)))
                    for a,b in [(0,2),(1,3)]]
    def free(v):return (u@v.reshape(16,16)@u.T).reshape(256)
    idx=np.arange(256);rbit=(idx>>3)&1;lbit=(idx>>7)&1
    before=free(vector);pair=np.zeros((2,256,256),complex);weights=[]
    for m in (0,1):
        branch=before*(rbit==m);weights.append(float(np.vdot(branch,branch).real))
        branch=free(branch)
        for b in (0,1):
            retained=branch*(lbit==b)
            jump=np.zeros(256,complex)
            targets=idx[lbit==b]
            jump[targets]=branch[targets^128]
            pair[b]+=np.outer(retained,retained.conj())+np.outer(jump,jump.conj())
    return pair[1]-pair[0],dict(thermal_marginal_errors=thermal_errors,early_weights=weights,
        pair_trace_error=float(np.max(abs(np.trace(pair,axis1=-2,axis2=-1)-1))),
        minimum_source_eigenvalue=float(np.linalg.eigvalsh(pair).min()))

def unitary(rho,u,sites,n):
    """Explicit local tensor contraction rather than bit-permutation updates."""
    outside=[q for q in range(n) if q not in sites]
    order=sites+outside+[n+q for q in sites+outside]
    ds,do=2**len(sites),2**len(outside)
    block=rho.reshape([2]*(2*n)).transpose(order).reshape(ds,do,ds,do)
    out=np.einsum('ai,ixjy,bj->axby',u,block,u.conj(),optimize=True)
    return out.reshape([2]*(2*n)).transpose(np.argsort(order)).reshape(2**n,2**n)

def relaxation(rho,us,n,resource=64):
    """Sequential local population/coherence updates from the CPTP equations."""
    population=math.exp(-us/(400*resource))
    coherence=math.exp(-us/(200*resource))
    out=rho.copy();view=out.reshape([2]*(2*n))
    for q in range(n):
        slices={}
        for a,b in itertools.product(range(2),repeat=2):
            key=[slice(None)]*(2*n);key[q]=a;key[n+q]=b
            slices[a,b]=tuple(key)
        excited=view[slices[1,1]].copy()
        view[slices[0,0]]+=(1-population)*excited
        view[slices[1,1]]=population*excited
        view[slices[0,1]]*=coherence
        view[slices[1,0]]*=coherence
    return out

def replacement(rho,sites,n):
    """Direct partial trace and tensor insertion of maximally mixed sites."""
    outside=[q for q in range(n) if q not in sites]
    order=outside+sites+[n+q for q in outside+sites]
    ds,do=2**len(sites),2**len(outside)
    block=rho.reshape([2]*(2*n)).transpose(order).reshape(do,ds,do,ds)
    reduced=np.trace(block,axis1=1,axis2=3)
    out=np.einsum('ab,ij->aibj',reduced,np.eye(ds)/ds)
    return out.reshape([2]*(2*n)).transpose(np.argsort(order)).reshape(2**n,2**n)

def gate_matrix(g,fraction=1.):
    if g['name']=='sx':
        angle=math.pi*fraction/4
        return math.cos(angle)*I-1j*math.sin(angle)*X
    if g['name']=='rz':
        angle=g['angle']*fraction/2
        return np.diag([np.exp(-1j*angle),np.exp(1j*angle)])
    if g['name']=='cz':return np.diag([1,1,1,np.exp(-1j*math.pi*fraction)])
    raise ValueError(g)

def interpreter_checks():
    rng=np.random.default_rng(240024)
    x=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho=x@x.conj().T;rho/=np.trace(rho);errors=[]
    for us in (.0175,.1,2.3):
        p=math.exp(-us/(400*64));deph=math.exp(-us/(800/3*64))
        amplitude=[np.diag([1,math.sqrt(p)]),np.array([[0,math.sqrt(1-p)],[0,0]])]
        phase=[math.sqrt((1+deph)/2)*I,math.sqrt((1-deph)/2)*Z]
        explicit=rho.copy()
        for q in range(3):
            explicit=sum(unitary(explicit,a@b,[q],3) for a in phase for b in amplitude)
        errors.append(float(np.max(abs(explicit-relaxation(rho,us,3)))))
    for sites in ([0],[2],[0,2]):
        ref=np.zeros((8,8),complex)
        ds=2**len(sites)
        for i,j in itertools.product(range(ds),repeat=2):
            k=np.zeros((ds,ds),complex);k[i,j]=1/math.sqrt(ds)
            ref+=unitary(rho,k,list(sites),3)
        errors.append(float(np.max(abs(ref-replacement(rho,list(sites),3)))))
    assert max(errors)<1e-12
    return max(errors)

def transport(source,timeline):
    events={}
    for g in timeline:
        events.setdefault(g['start_ns'],{'start':[],'end':[],'virtual':[]})
        if g['name']=='rz':events[g['start_ns']]['virtual'].append(g)
        else:
            events[g['start_ns']]['start'].append(g)
            events.setdefault(g['end_ns'],{'start':[],'end':[],'virtual':[]})['end'].append(g)
    times=sorted(events);active={};rho=source.copy();checkpoint=0
    for k,ns in enumerate(times):
        event=events[ns]
        for g in event['end']:
            assert g['index'] in active;del active[g['index']]
            p=(.0003 if g['name']=='sx' else .003)/64
            rho=(1-p)*rho+p*replacement(rho,g['qubits'],8)
        for g in event['virtual']:rho=unitary(rho,gate_matrix(g),g['qubits'],8)
        for g in event['start']:
            assert not any(set(g['qubits'])&set(other['qubits']) for other in active.values())
            active[g['index']]=g
        if k==len(times)-1:break
        dt=times[k+1]-ns
        rho=relaxation(rho,dt/2000,8)
        for g in active.values():
            fraction=dt/(g['end_ns']-g['start_ns'])
            rho=unitary(rho,gate_matrix(g,fraction),g['qubits'],8)
        rho=relaxation(rho,dt/2000,8)
        if times[k+1]-checkpoint>=100000:
            checkpoint=times[k+1]
            print('ALTERNATE full7001/64 completed_us',checkpoint/1000,flush=True)
    assert not active and times[-1]==229460
    return rho,dict(complete_transport=True,scheduled_us=times[-1]/1000,event_intervals=len(times)-1)

def assignment():
    return kron([np.array([[.99,.02],[.01,.98]])]*4).real

def receiver(delta,policy):
    rho=delta.copy();latch=0.
    for g in policy['gates']:
        if g['name']=='rz':rho=unitary(rho,gate_matrix(g),g['qubits'],4)
        else:
            us=.035 if g['name']=='sx' else .2
            rho=relaxation(rho,us/2,4)
            rho=unitary(rho,gate_matrix(g),g['qubits'],4)
            rho=relaxation(rho,us/2,4)
            p=(.0003 if g['name']=='sx' else .003)/64
            rho=(1-p)*rho+p*replacement(rho,g['qubits'],4);latch+=us
    rho=relaxation(rho,.5,4);latch+=.5
    reported=assignment()@np.diag(rho).real
    score=.5+.5*reported[np.array(policy['decoder'])==1].sum()
    return float(score),rho,reported,latch

def estimator(counts):
    true=np.linalg.solve(assignment(),(counts/2048).reshape(-1,16).T).T.reshape(2,81,16)
    settings=list(itertools.product('XYZ',repeat=4));words=list(itertools.product('IXYZ',repeat=4))
    labels=((np.arange(16)[:,None]>>np.arange(3,-1,-1))&1)
    pauli={'I':I,'X':X,'Y':Y,'Z':Z};linear=np.zeros((2,16,16),complex)
    for word in words:
        compatible=[k for k,setting in enumerate(settings) if all(a=='I' or a==b for a,b in zip(word,setting))]
        occupied=[q for q,a in enumerate(word) if a!='I']
        signs=(-1.)**labels[:,occupied].sum(axis=1)
        coefficient=(true[:,compatible,:]@signs).mean(axis=1)
        if not occupied:coefficient=np.ones(2)
        linear+=coefficient[:,None,None]*kron([pauli[a] for a in word])/16
    output=[]
    for rho in linear:
        eigenvalues,basis=np.linalg.eigh((rho+rho.conj().T)/2)
        ordered=np.sort(eigenvalues)[::-1]
        for k in range(1,17):
            threshold=(ordered[:k].sum()-1)/k
            if k==16 or ordered[k]<=threshold:break
        weights=np.maximum(eigenvalues-threshold,0)
        output.append((basis*weights)@basis.conj().T)
    return np.array(output)

def load_npz(name):
    content=base64.b64decode((DATA/(name+'.b64')).read_text(),validate=False)
    return np.load(io.BytesIO(content),allow_pickle=False)

def main():
    started=time.perf_counter()
    initial,source_status=source(7001)
    with load_npz('source_seed_7001.npz') as z:source_error=float(np.max(abs(initial-z['difference'])))
    assert source_error<1e-9 and max(source_status['thermal_marginal_errors'])<1e-9
    assert source_status['pair_trace_error']<1e-9 and source_status['minimum_source_eigenvalue']>-1e-9
    schedule=base64.b64decode((DATA/'schedule_seed_7001.json.gz.b64').read_text())
    timeline=json.loads(gzip.decompress(schedule))
    full,status=transport(initial,timeline)
    right=np.trace(full.reshape(16,16,16,16),axis1=0,axis2=2)
    with load_npz('physical64_seed_7001.npz') as z:state_error=float(np.max(abs(right-z['difference'])))
    assert state_error<1e-9
    references=json.loads((DATA/'results64_seed_7001.json').read_text())
    cases=[]
    for label in ['oracle','n2048_d0_nominal','n2048_d1_nominal','n2048_d2_nominal']:
        policy=json.loads((DATA/f'policy_{label}_seed_7001.json').read_text())
        score,rho,reported,latch=receiver(right,policy)
        reference=next(x for x in references if x['label']==label)
        with load_npz(f'result64_{label}_seed_7001.npz') as z:
            receiver_error=float(np.max(abs(rho-z['receiver_difference'])))
            readout_error=float(np.max(abs(reported-z['reported_difference'])))
        assert max(abs(score-reference['success']),receiver_error,readout_error)<1e-9
        estimator_error=None
        if label!='oracle':
            with load_npz(f'estimate_{label}_seed_7001.npz') as z:
                counts=z['counts'];pair=z['estimated_pair']
            estimator_error=float(np.max(abs(estimator(counts)-pair)))
            assert estimator_error<1e-9
            with load_npz('tomography_seed_7001.npz') as z:probs=z['probabilities']
            draw=int(label.split('_')[1][1:])
            rng=np.random.default_rng(210021+100*7001+draw+10000000*2048)
            resampled=np.array([[rng.multinomial(2048,p) for p in bit] for bit in probs])
            assert np.array_equal(counts,resampled)
        cases.append({'label':label,'expected_success_probability':score,'reference_error':abs(score-reference['success']),
                      'receiver_error':receiver_error,'readout_error':readout_error,'estimator_error':estimator_error,'latch_us':latch})
    result={'status':'PASS_FULL_7001_FROZEN_CASE','source_error':source_error,'right_state_error':state_error,
            'transport':status,'cases':cases,'seconds':time.perf_counter()-started,
            'scope':'Same-project alternate equations on a frozen schedule/receiver, not independent schedule synthesis or external review.',
            'hardware_executions':0,'external_review_established':False,'numpy_version':np.__version__}
    output=HERE/'results';output.mkdir(exist_ok=True)
    (output/'full_case.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
