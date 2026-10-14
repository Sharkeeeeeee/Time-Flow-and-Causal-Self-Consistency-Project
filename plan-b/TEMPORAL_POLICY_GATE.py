"""Declared-time preflight only. Run alongside packet/hash/count validation.

No provider authentication, stationary-noise certificate or hardware submission.
The valid window and maximum age are registered operational assumptions and
must not be represented as measured device calibration lifetimes.
"""
from pathlib import Path
from datetime import datetime, timedelta
import json
import math

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def utc(s):
    d=datetime.fromisoformat(s.replace('Z','+00:00'))
    if d.tzinfo is None or d.utcoffset()!=timedelta(0):raise ValueError('UTC required')
    return d

def temporal_gate(binding_path,counts_path,policy):
    directory=Path(binding_path).parent
    snapshot=read(directory/'snapshot.json');certificate=read(directory/'calibration_handoff.json')
    rows=read(counts_path)['rows']
    if not rows:raise ValueError('No execution rows')
    if policy.get('calibration_id')!=certificate['calibration_id']:raise ValueError('Policy calibration mismatch')
    lower,upper=utc(policy['calibration_valid_from_utc']),utc(policy['calibration_valid_until_utc'])
    age=policy['snapshot_max_age_seconds']
    if isinstance(age,bool) or not isinstance(age,(int,float)) or not math.isfinite(age) or age<=0:
        raise ValueError('Finite positive maximum snapshot age required')
    if upper<lower:raise ValueError('Invalid calibration window')
    fetched=[utc(snapshot['fetched_utc']),utc(certificate['fetched_utc'])]
    for row in rows:
        start,end=utc(row['executed_start_utc']),utc(row['executed_end_utc'])
        if end<start:raise ValueError('Execution end precedes start')
        if start<lower or end>upper:raise ValueError('Execution outside registered calibration window')
        if any(start<t for t in fetched):raise ValueError('Execution predates supporting snapshot/handoff')
        if any((end-t).total_seconds()>age for t in fetched):raise ValueError('Snapshot/handoff too old')
    return {'status':'PASS_DECLARED_TIME_CONSISTENCY_ONLY','rows':len(rows)}

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['binding','counts','policy']:parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(temporal_gate(args.binding,args.counts,read(args.policy)),indent=2))
