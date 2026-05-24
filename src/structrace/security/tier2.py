from __future__ import annotations

import base64
import csv
import getpass
import hashlib
import platform
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Union

from cryptography.fernet import Fernet, InvalidToken

from .common import sha256_file


def machine_fingerprint() -> str:
    raw = f"{platform.node()}|{platform.system()}|{platform.machine()}|{uuid.getnode()}|{getpass.getuser()}"
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _fernet_from_secret(secret: str) -> Fernet:
    digest = hashlib.sha256(secret.encode("utf-8")).digest()
    return Fernet(base64.urlsafe_b64encode(digest))


def encrypt_file(input_path: Union[str, Path], output_path: Union[str, Path], *, machine_id: Optional[str] = None) -> str:
    machine_id = machine_id or machine_fingerprint()
    token = _fernet_from_secret(machine_id).encrypt(Path(input_path).read_bytes())
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_bytes(token)
    return machine_id


def decrypt_file(input_path: Union[str, Path], output_path: Union[str, Path], *, machine_id: Optional[str] = None) -> None:
    machine_id = machine_id or machine_fingerprint()
    try:
        data = _fernet_from_secret(machine_id).decrypt(Path(input_path).read_bytes())
    except InvalidToken as exc:
        raise ValueError("Decryption failed: machine_id does not match or ciphertext is invalid.") from exc
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_bytes(data)


def audit_artifacts(paths: List[Union[str, Path]], output_csv: Optional[Union[str, Path]] = None) -> List[Dict[str, Union[str, int]]]:
    rows = [
        {
            "artifact": str(path),
            "bytes": Path(path).stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
    ]
    if output_csv is not None:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        with Path(output_csv).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    return rows
