from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional

from .robustness import add_gaussian_noise_to_atoms, round_pdb_coordinates, translate_pdb
from .security import (
    LedgerEvent,
    audit_artifacts,
    build_ledger,
    build_public_registry,
    decrypt_file,
    encrypt_file,
    machine_fingerprint,
    verify_ledger,
    verify_public_registry,
)
from .watermark import decode_bits, decode_text, embed_bits, embed_text
from .watermark.payload import text_to_bits


def _write_csv(path: str, rows: List[Dict]) -> None:
    if not rows:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="structrace", description="StrucTrace watermarking and safeguard CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    embed = sub.add_parser("embed", help="Embed a text payload or bit string into a PDB file")
    embed.add_argument("input_pdb")
    embed.add_argument("-o", "--output", required=True)
    embed.add_argument("--text")
    embed.add_argument("--bits")
    embed.add_argument("--start-freq", type=int, default=5)
    embed.add_argument("--strength", type=float, default=0.02)

    decode = sub.add_parser("decode", help="Decode a watermark from a query PDB using a master PDB")
    decode.add_argument("master_pdb")
    decode.add_argument("query_pdb")
    decode.add_argument("--bits", type=int, required=True)
    decode.add_argument("--expected-text")
    decode.add_argument("--expected-bits")
    decode.add_argument("--start-freq", type=int, default=5)

    attack = sub.add_parser("attack", help="Generate simple coordinate perturbations")
    attack_sub = attack.add_subparsers(dest="attack_command", required=True)
    translate = attack_sub.add_parser("translate")
    translate.add_argument("input_pdb")
    translate.add_argument("-o", "--output", required=True)
    translate.add_argument("--vector", nargs=3, type=float, required=True)
    rounding = attack_sub.add_parser("round")
    rounding.add_argument("input_pdb")
    rounding.add_argument("-o", "--output", required=True)
    rounding.add_argument("--decimals", type=int, default=3)
    noise = attack_sub.add_parser("noise")
    noise.add_argument("input_pdb")
    noise.add_argument("-o", "--output", required=True)
    noise.add_argument("--sigma", type=float, required=True)
    noise.add_argument("--seed", type=int, default=1)

    tier1 = sub.add_parser("tier1", help="Public provenance registry utilities")
    tier1_sub = tier1.add_subparsers(dest="tier1_command", required=True)
    registry = tier1_sub.add_parser("build-registry")
    registry.add_argument("--record", action="append", nargs=4, metavar=("CASE_ID", "MASTER_PDB", "QUERY_PDB", "PAYLOAD"), required=True)
    registry.add_argument("-o", "--output", required=True)
    verify_registry = tier1_sub.add_parser("verify")
    verify_registry.add_argument("registry_csv")
    verify_registry.add_argument("-o", "--output")

    tier2 = sub.add_parser("tier2", help="Hardware-bound encryption and artifact audit")
    tier2_sub = tier2.add_subparsers(dest="tier2_command", required=True)
    fp = tier2_sub.add_parser("fingerprint")
    fp.set_defaults(_fingerprint=True)
    encrypt = tier2_sub.add_parser("encrypt")
    encrypt.add_argument("input_file")
    encrypt.add_argument("-o", "--output", required=True)
    encrypt.add_argument("--machine-id")
    decrypt = tier2_sub.add_parser("decrypt")
    decrypt.add_argument("input_file")
    decrypt.add_argument("-o", "--output", required=True)
    decrypt.add_argument("--machine-id")
    audit = tier2_sub.add_parser("audit")
    audit.add_argument("paths", nargs="+")
    audit.add_argument("-o", "--output")

    tier3 = sub.add_parser("tier3", help="DRM ledger utilities")
    tier3_sub = tier3.add_subparsers(dest="tier3_command", required=True)
    ledger = tier3_sub.add_parser("build-ledger")
    ledger.add_argument("--event", action="append", nargs=5, metavar=("ASSET_ID", "EVENT", "ACTOR", "COUNTERPARTY", "PERMISSION"), required=True)
    ledger.add_argument("--json", required=True)
    ledger.add_argument("--csv")
    ledger_verify = tier3_sub.add_parser("verify")
    ledger_verify.add_argument("ledger_json")

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "embed":
        if bool(args.text) == bool(args.bits):
            raise SystemExit("Provide exactly one of --text or --bits.")
        result = (
            embed_text(args.input_pdb, args.text, args.output, start_freq=args.start_freq, strength=args.strength)
            if args.text is not None
            else embed_bits(args.input_pdb, args.bits, args.output, start_freq=args.start_freq, strength=args.strength)
        )
        print(json.dumps(result.__dict__, indent=2))
        return

    if args.command == "decode":
        if args.expected_text is not None:
            result = decode_text(args.master_pdb, args.query_pdb, args.bits, expected_text=args.expected_text, start_freq=args.start_freq)
        else:
            result = decode_bits(args.master_pdb, args.query_pdb, args.bits, expected_bits=args.expected_bits, start_freq=args.start_freq)
        print(json.dumps(result.__dict__, indent=2))
        return

    if args.command == "attack":
        if args.attack_command == "translate":
            translate_pdb(args.input_pdb, args.output, tuple(args.vector))
        elif args.attack_command == "round":
            round_pdb_coordinates(args.input_pdb, args.output, args.decimals)
        elif args.attack_command == "noise":
            add_gaussian_noise_to_atoms(args.input_pdb, args.output, args.sigma, seed=args.seed)
        return

    if args.command == "tier1":
        if args.tier1_command == "build-registry":
            records = [(case, master, query, payload) for case, master, query, payload in args.record]
            rows = build_public_registry(records, args.output)
            print(json.dumps([row.__dict__ for row in rows], indent=2))
        elif args.tier1_command == "verify":
            rows = verify_public_registry(args.registry_csv)
            if args.output:
                _write_csv(args.output, rows)
            print(json.dumps(rows, indent=2))
        return

    if args.command == "tier2":
        if args.tier2_command == "fingerprint":
            print(machine_fingerprint())
        elif args.tier2_command == "encrypt":
            machine_id = encrypt_file(args.input_file, args.output, machine_id=args.machine_id)
            print(json.dumps({"output": args.output, "machine_id": machine_id}, indent=2))
        elif args.tier2_command == "decrypt":
            decrypt_file(args.input_file, args.output, machine_id=args.machine_id)
            print(json.dumps({"output": args.output}, indent=2))
        elif args.tier2_command == "audit":
            rows = audit_artifacts(args.paths, args.output)
            print(json.dumps(rows, indent=2))
        return

    if args.command == "tier3":
        if args.tier3_command == "build-ledger":
            events = [LedgerEvent(asset_id, event, actor, counterparty, permission) for asset_id, event, actor, counterparty, permission in args.event]
            rows = build_ledger(events, output_json=args.json, output_csv=args.csv)
            print(json.dumps(rows, indent=2))
        elif args.tier3_command == "verify":
            print(json.dumps({"valid": verify_ledger(args.ledger_json)}, indent=2))
        return


if __name__ == "__main__":
    main()
