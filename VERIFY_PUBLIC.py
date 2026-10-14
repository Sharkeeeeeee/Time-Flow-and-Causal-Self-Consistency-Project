"""SPDX-License-Identifier: MIT. Offline integrity only, not external scientific review."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parent
manifest=json.loads((root/'PUBLIC_MANIFEST.json').read_text())
failures=[]
for name,expected in manifest.items():
    path=(root/name).resolve()
    if not path.is_relative_to(root) or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
        failures.append(name)
print(json.dumps({'status':'FAIL' if failures else 'PASS_INTEGRITY','files_checked':len(manifest),'failures':failures,
                  'hardware_executions':0,'external_review_established':False},indent=2))
raise SystemExit(bool(failures))
