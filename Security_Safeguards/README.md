# Security safeguards

This directory provides compact scripts and case outputs for the three-tier safeguard model shown in Figure 5.

## Structure

- `Tier1_Public_Provenance/`: public registry and verification example. This tier records a structure fingerprint and provenance payload so that a submitted structure can be matched to a registered master record before reference-guided decoding.
- `Tier2_Hardware_Bound_Access/`: hardware-bound access case. This tier preserves the 1BIL encrypted/decrypted example and audits file integrity for authorized versus protected structural assets.
- `Tier3_Digital_Rights_Management/`: digital-rights ledger example. This tier records a minimal chain of licensing transactions for auditable ownership, permission and revocation events.

The scripts are intentionally lightweight and use only the Python standard library. They are not production security services; they are reproducible evidence of the workflow logic used to support the three safeguards.
