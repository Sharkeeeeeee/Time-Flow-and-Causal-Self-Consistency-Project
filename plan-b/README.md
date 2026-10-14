# Plan B: bounded diagnostics evidence

The existing utilities cover planning arithmetic and a declared-time consistency gate. They do not constitute the proposed standalone analyzer, two non-SYK user examples, power study, interoperability or external adoption.

```text
python -I -S -B plan-b/REPRODUCE_COSTS.py
python -I -S -B plan-b/TEST_TEMPORAL_POLICY.py
```

Only Python's standard library is required. The first command shows why the historical conservative four-outcome drift rule cannot pass a 0.10 tolerance at 4,096 shots/task: its paired uncertainty radius is about 0.1681. The necessary width threshold at zero observed discrepancy is 11,575 shots, which is not a power guarantee. The planning prices use published rates checked on 14 October 2026; they are not quotes, purchases or observed savings.

The second command accepts one valid synthetic case and rejects eight invalid cases, including timestamps before supporting snapshots, expired windows and stale supporting records. A declared window is an operational policy, not a measured device lifetime. Hashes and timestamps alone do not authenticate a provider or establish stationarity.

Both commands run offline. A third party can retain their stdout, environment and findings, but a successful automated run does not establish reviewer identity or independence. Hardware executions and actual external review responses remain zero.

Price reference: [AWS Braket](https://aws.amazon.com/braket/pricing/). Statistical and physical claims must remain separate from dollar amounts.
