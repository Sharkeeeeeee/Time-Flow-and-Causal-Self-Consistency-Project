"""Synthetic standard-library regression of the declared-time preflight."""
from pathlib import Path
import importlib.util
import json
import tempfile

def main():
    base=Path(__file__).resolve().parent
    core=base/'TEMPORAL_POLICY_GATE.py'
    if not core.exists():core=base/'pilot_temporal_policy_v31.py'
    spec=importlib.util.spec_from_file_location('time_gate',core)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    policy={'calibration_id':'SYNTHETIC','calibration_valid_from_utc':'2026-10-14T00:00:00Z',
            'calibration_valid_until_utc':'2026-10-14T01:00:00Z','snapshot_max_age_seconds':3600}
    cases=[('valid','2026-10-14T00:01:00Z','2026-10-14T00:01:01Z',{},True),
      ('before_snapshot','2026-10-13T00:01:00Z','2026-10-13T00:01:01Z',{'calibration_valid_from_utc':'2026-10-01T00:00:00Z'},False),
      ('outside_window','2026-10-14T02:00:00Z','2026-10-14T02:00:01Z',{},False),
      ('old_snapshot','2026-10-14T02:00:00Z','2026-10-14T02:00:01Z',{'calibration_valid_until_utc':'2026-10-15T00:00:00Z'},False),
      ('expires_during_run','2026-10-14T00:59:59Z','2026-10-14T01:00:01Z',{'calibration_valid_until_utc':'2026-10-15T00:00:00Z'},False),
      ('reversed_times','2026-10-14T00:01:01Z','2026-10-14T00:01:00Z',{},False),
      ('missing_utc','2026-10-14T00:01:00','2026-10-14T00:01:01',{},False),
      ('bad_age','2026-10-14T00:01:00Z','2026-10-14T00:01:01Z',{'snapshot_max_age_seconds':float('nan')},False),
      ('wrong_calibration','2026-10-14T00:01:00Z','2026-10-14T00:01:01Z',{'calibration_id':'different'},False)]
    output={}
    with tempfile.TemporaryDirectory(prefix='temporal_gate_v31_',dir=base) as tmp:
        d=Path(tmp)
        assert d.resolve().parent==base.resolve()
        def write(name,value):(d/name).write_text(json.dumps(value)+'\n',encoding='utf-8')
        write('snapshot.json',{'fetched_utc':'2026-10-14T00:00:00Z'})
        write('calibration_handoff.json',{'calibration_id':'SYNTHETIC','fetched_utc':'2026-10-14T00:00:00Z'})
        for name,start,end,patch,expected in cases:
            write('counts.json',{'rows':[{'executed_start_utc':start,'executed_end_utc':end}]})
            try:module.temporal_gate(d/'binding.json',d/'counts.json',{**policy,**patch});passed=True
            except (ValueError,KeyError):passed=False
            assert passed==expected,name
            output[name]='ACCEPTED' if passed else 'REJECTED'
    print(json.dumps({'status':'PASS_SYNTHETIC_TEMPORAL_TESTS','cases':output,'hardware_executions':0},indent=2))

if __name__=='__main__':main()
