#!/usr/bin/env python
"""
Upload files to an existing Zenodo draft, safely re-runnable.

  read -s -p "Zenodo token: " ZENODO_TOKEN; export ZENODO_TOKEN; echo
  python zenodo_upload.py --deposition 1234567 final/*.nc final/MANIFEST.txt README.md

For each file:
  - already in the draft with the same size (and MD5, unless --no-md5) -> skipped
  - otherwise uploaded (replacing any partial/mismatched copy), and the MD5
    Zenodo reports is checked against the local file.
So after an interruption, just run the same command again.

--sandbox targets sandbox.zenodo.org (separate account and token) for a
rehearsal. The script never publishes: review the draft in the browser and
press Publish there.
"""
import argparse
import hashlib
import os
import sys
import time

import requests


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def human(n):
    return f"{n / 1e9:.2f} GB" if n >= 1e8 else f"{n / 1e6:.1f} MB"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--deposition", required=True, help="draft id (number in the upload URL)")
    ap.add_argument("--sandbox", action="store_true")
    ap.add_argument("--no-md5", action="store_true",
                    help="skip local checksums (faster; compares sizes only)")
    a = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("set ZENODO_TOKEN first (see the usage at the top of this script)")
    host = "https://sandbox.zenodo.org" if a.sandbox else "https://zenodo.org"
    auth = {"Authorization": f"Bearer {token}"}

    r = requests.get(f"{host}/api/deposit/depositions/{a.deposition}", headers=auth, timeout=60)
    if r.status_code != 200:
        sys.exit(f"could not read draft {a.deposition}: HTTP {r.status_code} {r.text[:300]}")
    dep = r.json()
    if dep.get("submitted"):
        sys.exit("this record is already published; create a new version for new files")
    bucket = dep["links"]["bucket"]
    remote = {f["filename"]: (int(f["filesize"]), f["checksum"].replace("md5:", ""))
              for f in dep.get("files", [])}

    missing = [f for f in a.files if not os.path.isfile(f)]
    if missing:
        sys.exit("not found: " + ", ".join(missing))
    total = sum(os.path.getsize(f) for f in a.files)
    print(f"draft {a.deposition} on {host}: {len(remote)} file(s) already there; "
          f"{len(a.files)} local file(s), {human(total)}\n")

    failed = []
    for path in a.files:
        name, size = os.path.basename(path), os.path.getsize(path)
        local = None if a.no_md5 else md5(path)
        if name in remote and remote[name][0] == size and (a.no_md5 or remote[name][1] == local):
            print(f"skip      {name} (already uploaded)")
            continue
        print(f"upload    {name} ({human(size)}) ...", end="", flush=True)
        t0 = time.time()
        try:
            with open(path, "rb") as fh:
                r = requests.put(f"{bucket}/{name}", data=fh, headers=auth, timeout=None)
            r.raise_for_status()
            got = r.json().get("checksum", "").replace("md5:", "")
        except Exception as e:
            print(f" FAILED: {e}")
            failed.append(name)
            continue
        rate = size / max(time.time() - t0, 1e-3) / 1e6
        if local and got != local:
            print(f" CHECKSUM MISMATCH (local {local}, Zenodo {got})")
            failed.append(name)
        else:
            print(f" ok, {rate:.0f} MB/s" + ("" if a.no_md5 else ", md5 verified"))

    print()
    if failed:
        print(f"{len(failed)} file(s) failed: {', '.join(failed)}\nRe-run the same command to retry.")
        sys.exit(1)
    print(f"All files are in the draft. Review and publish at {host}/uploads/{a.deposition}")


if __name__ == "__main__":
    main()
