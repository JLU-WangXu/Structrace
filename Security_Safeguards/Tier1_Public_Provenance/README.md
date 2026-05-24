# Tier 1: Public provenance

This tier demonstrates a public provenance registry for reference-guided verification. The registry stores a master structure fingerprint, a watermarked query fingerprint, the selected C-alpha table and the decoded provenance payload.

Run:

```bash
python public_provenance_verification.py
```

Outputs:

- `public_registry.csv`: compact registry records for 6MRR, 8HFE and 8VC8.
- `public_verification_results.csv`: integrity check results for registered watermarked structures.
