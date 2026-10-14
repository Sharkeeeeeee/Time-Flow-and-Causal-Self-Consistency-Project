from pathlib import Path
import csv,json,statistics
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'data/results.csv').open(newline='')))
selected=[r for r in rows if r['resource']=='64' and r['label'] in ['n2048_d0_nominal','n2048_d1_nominal','n2048_d2_nominal']]
assert len(rows)==408 and len(selected)==36
groups={seed:[float(r['success']) for r in selected if r['seed']==seed] for seed in sorted({r['seed'] for r in selected})}
assert len(groups)==12 and all(len(x)==3 for x in groups.values())
means=[statistics.mean(x) for x in groups.values()]
result={'rows':len(rows),'models':len(groups),'draws_per_model':3,'mean':statistics.mean(means),
        'minimum_model_mean':min(means),'minimum_draw':min(float(r['success']) for r in selected),
        'model_means_at_least_70_percent':sum(x>=.70 for x in means),
        'scope':'Arithmetic from archived simulator probabilities; not a fresh 12-model forward simulation.'}
assert abs(result['mean']-.712627382528563)<1e-12
assert abs(result['minimum_model_mean']-.6942280612527658)<1e-12
assert result['model_means_at_least_70_percent']==8
print(json.dumps(result,indent=2))
