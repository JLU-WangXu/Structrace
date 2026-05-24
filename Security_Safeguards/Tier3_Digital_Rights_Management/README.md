# Tier 3: Digital rights management

This tier provides a minimal auditable ledger for licensing events. Each event is chained to the previous record hash, illustrating how structural-asset permissions can be tracked for licensing, collaboration and revocation.

Run:

```bash
python drm_ledger_demo.py
```

Outputs:

- `drm_ledger.json`: structured chained ledger records.
- `drm_ledger.csv`: tabular view of the same events.
