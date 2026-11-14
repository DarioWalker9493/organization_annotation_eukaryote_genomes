#!/usr/bin/env python3

import re
from pathlib import Path

FINALDIR = Path("output/MAKER/final")

PROTEIN_FASTA = FINALDIR / "assembly.all.maker.proteins.fasta.renamed.fasta"
IPRSCAN = FINALDIR / "output.iprscan"

AED_THRESHOLD = 0.5


def main():
    if not PROTEIN_FASTA.exists():
        raise FileNotFoundError(PROTEIN_FASTA)

    if not IPRSCAN.exists():
        raise FileNotFoundError(IPRSCAN)

    # Proteins with at least one Pfam hit
    pfam_ids = set()

    with IPRSCAN.open() as handle:
        for line in handle:
            if line.strip():
                pfam_ids.add(line.split("\t", 1)[0])

    stats = {
        "low": {"total": 0, "pfam": 0},
        "high": {"total": 0, "pfam": 0},
    }

    total_proteins = 0
    proteins_with_aed = 0

    with PROTEIN_FASTA.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            total_proteins += 1

            protein_id = line[1:].split()[0]

            match = re.search(r"\bAED:([0-9.]+)", line)
            if match is None:
                continue

            proteins_with_aed += 1
            aed = float(match.group(1))

            group = "low" if aed <= AED_THRESHOLD else "high"

            stats[group]["total"] += 1

            if protein_id in pfam_ids:
                stats[group]["pfam"] += 1

    print("========================================")
    print("Week 3 - AED / Pfam analysis")
    print("========================================")
    print()
    print(f"Total protein records:          {total_proteins}")
    print(f"Proteins with AED value:        {proteins_with_aed}")
    print(f"Proteins with >=1 Pfam hit:     {len(pfam_ids)}")

    if total_proteins:
        overall_pct = 100 * len(pfam_ids) / total_proteins
        print(f"Overall Pfam coverage:          {overall_pct:.2f}%")

    print()
    print(f"AED threshold: {AED_THRESHOLD:.2f}")
    print()

    for label, description in [
        ("low", f"AED <= {AED_THRESHOLD:.2f}"),
        ("high", f"AED > {AED_THRESHOLD:.2f}"),
    ]:
        total = stats[label]["total"]
        pfam = stats[label]["pfam"]
        pct = 100 * pfam / total if total else 0

        print(
            f"{description}: "
            f"{pfam}/{total} proteins with Pfam = {pct:.2f}%"
        )

    print()
    print("========================================")


if __name__ == "__main__":
    main()
