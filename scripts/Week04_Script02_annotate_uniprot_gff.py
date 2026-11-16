#!/usr/bin/env python3

import argparse
import re
from urllib.parse import quote


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add UniProt BLAST annotations to MAKER GFF3 and summarize hit/no-hit proteins."
    )
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--besthits", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output-gff", required=True)
    parser.add_argument("--output-table", required=True)
    return parser.parse_args()


def parse_uniprot_fasta(path):
    lookup = {}

    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            parts = header.split(None, 1)

            token = parts[0]
            description = parts[1] if len(parts) > 1 else ""

            # Keep protein function only, remove organism metadata.
            if " OS=" in description:
                description = description.split(" OS=", 1)[0].strip()

            keys = {token}

            # Standard UniProt header:
            # sp|P12345|ENTRY_NAME
            token_parts = token.split("|")

            if len(token_parts) >= 3:
                keys.add(token_parts[1])
                keys.add(token_parts[2])

            # Also support possible accession versions.
            for key in list(keys):
                keys.add(re.sub(r"\.\d+$", "", key))

            record = {
                "token": token,
                "accession": token_parts[1] if len(token_parts) >= 3 else token,
                "function": description
            }

            for key in keys:
                lookup[key] = record

    return lookup


def parse_proteins(path):
    proteins = {}
    current = None

    with open(path) as handle:
        for line in handle:
            line = line.rstrip()

            if line.startswith(">"):
                current = line[1:].split()[0]
                proteins[current] = 0
            elif current:
                proteins[current] += len(line.strip())

    return proteins


def resolve_uniprot(subject, lookup):
    candidates = [subject]

    candidates.append(re.sub(r"\.\d+$", "", subject))

    parts = subject.split("|")

    if len(parts) >= 3:
        candidates.extend([
            parts[1],
            parts[2],
            re.sub(r"\.\d+$", "", parts[1]),
            re.sub(r"\.\d+$", "", parts[2])
        ])

    for key in candidates:
        if key in lookup:
            return lookup[key]

    return None


def classify_function(function):
    text = function.lower()

    uncharacterized_terms = [
        "uncharacterized",
        "hypothetical protein",
        "unknown protein",
        "protein of unknown function",
        "unnamed protein"
    ]

    for term in uncharacterized_terms:
        if term in text:
            return "uncharacterized"

    return "characterized"


def parse_besthits(path, uniprot_lookup):
    hits = {}

    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 12:
                continue

            query = fields[0]
            subject = fields[1]

            # File should already contain one best hit/query.
            if query in hits:
                continue

            record = resolve_uniprot(subject, uniprot_lookup)

            if record is None:
                hits[query] = {
                    "subject": subject,
                    "accession": "",
                    "function": "",
                    "function_class": "unresolved",
                    "evalue": fields[10],
                    "bitscore": fields[11]
                }
            else:
                hits[query] = {
                    "subject": subject,
                    "accession": record["accession"],
                    "function": record["function"],
                    "function_class": classify_function(record["function"]),
                    "evalue": fields[10],
                    "bitscore": fields[11]
                }

    return hits


def parse_attributes(text):
    attrs = {}

    for item in text.strip().strip(";").split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value

    return attrs


def encode_gff_value(value):
    return quote(
        value,
        safe="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._:-|+()/ "
    ).replace(" ", "%20")


def add_annotation_to_attributes(attr_text, hit):
    text = attr_text.rstrip().rstrip(";")

    if hit is None or not hit.get("function"):
        return text + ";"

    function = encode_gff_value(hit["function"])
    accession = hit["accession"]

    # Add UniProt cross-reference while retaining existing Dbxref.
    attrs = parse_attributes(text)

    if accession:
        uniprot_ref = "UniProtKB:" + accession

        if "Dbxref" in attrs:
            existing = attrs["Dbxref"].split(",")

            if uniprot_ref not in existing:
                old = "Dbxref=" + attrs["Dbxref"]
                new = old + "," + uniprot_ref
                text = text.replace(old, new, 1)
        else:
            text += ";Dbxref=" + uniprot_ref

    # maker_functional_gff normally stores functional information in Note.
    if "Note=" not in text:
        text += ";Note=" + function

    return text + ";"


def median(values):
    values = sorted(values)

    if not values:
        return 0

    n = len(values)
    middle = n // 2

    if n % 2:
        return values[middle]

    return (values[middle - 1] + values[middle]) / 2.0


def main():
    args = parse_args()

    uniprot = parse_uniprot_fasta(args.uniprot)
    proteins = parse_proteins(args.proteins)
    hits = parse_besthits(args.besthits, uniprot)

    # ------------------------------------------------------------------------
    # Parse Pfam presence from current final GFF.
    # ------------------------------------------------------------------------

    pfam = {}

    with open(args.gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            if fields[2] != "mRNA":
                continue

            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("ID")

            if transcript_id:
                pfam[transcript_id] = "Pfam:" in fields[8]

    # ------------------------------------------------------------------------
    # Write protein-level summary table.
    # ------------------------------------------------------------------------

    hit_lengths = []
    nohit_lengths = []

    hit_pfam = 0
    nohit_pfam = 0

    characterized = 0
    uncharacterized = 0
    unresolved = 0

    with open(args.output_table, "w") as out:
        out.write(
            "query_id\tlength_aa\thit_status\tuniprot_subject\t"
            "uniprot_accession\tfunction\tfunction_class\t"
            "evalue\tbitscore\tpfam_present\n"
        )

        for query in sorted(proteins):
            length = proteins[query]
            has_pfam = pfam.get(query, False)

            if query in hits:
                hit = hits[query]

                hit_lengths.append(length)

                if has_pfam:
                    hit_pfam += 1

                if hit["function_class"] == "characterized":
                    characterized += 1
                elif hit["function_class"] == "uncharacterized":
                    uncharacterized += 1
                else:
                    unresolved += 1

                values = [
                    query,
                    str(length),
                    "hit",
                    hit["subject"],
                    hit["accession"],
                    hit["function"].replace("\t", " "),
                    hit["function_class"],
                    hit["evalue"],
                    hit["bitscore"],
                    "yes" if has_pfam else "no"
                ]

            else:
                nohit_lengths.append(length)

                if has_pfam:
                    nohit_pfam += 1

                values = [
                    query,
                    str(length),
                    "no_hit",
                    "",
                    "",
                    "",
                    "no_hit",
                    "",
                    "",
                    "yes" if has_pfam else "no"
                ]

            out.write("\t".join(values) + "\n")

    # ------------------------------------------------------------------------
    # Prepare gene-level mapping.
    # Pick highest-bit-score transcript hit if a gene has multiple isoforms.
    # ------------------------------------------------------------------------

    gene_hits = {}

    for transcript, hit in hits.items():
        gene = re.sub(r"-R[A-Z]+$", "", transcript)

        try:
            bitscore = float(hit["bitscore"])
        except ValueError:
            bitscore = 0.0

        if gene not in gene_hits:
            gene_hits[gene] = (bitscore, hit)
        elif bitscore > gene_hits[gene][0]:
            gene_hits[gene] = (bitscore, hit)

    # ------------------------------------------------------------------------
    # Annotate final GFF3.
    # ------------------------------------------------------------------------

    annotated_genes = 0
    annotated_mrnas = 0

    with open(args.gff) as inp, open(args.output_gff, "w") as out:

        for line in inp:

            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                out.write(line)
                continue

            feature = fields[2]
            attrs = parse_attributes(fields[8])

            hit = None

            if feature == "mRNA":
                feature_id = attrs.get("ID")

                if feature_id in hits:
                    hit = hits[feature_id]
                    annotated_mrnas += 1

            elif feature == "gene":
                feature_id = attrs.get("ID")

                if feature_id in gene_hits:
                    hit = gene_hits[feature_id][1]
                    annotated_genes += 1

            if hit is not None:
                fields[8] = add_annotation_to_attributes(
                    fields[8],
                    hit
                )

            out.write("\t".join(fields) + "\n")

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    total = len(proteins)
    n_hits = len(hit_lengths)
    n_nohits = len(nohit_lengths)

    print("========================================")
    print("Week 4 UniProt functional annotation")
    print("========================================")
    print("")
    print("Total proteins:              {}".format(total))
    print("Proteins with UniProt hit:   {} ({:.2f}%)".format(
        n_hits,
        100.0 * n_hits / total
    ))
    print("Proteins without hit:        {} ({:.2f}%)".format(
        n_nohits,
        100.0 * n_nohits / total
    ))
    print("")
    print("Characterized UniProt hits:  {}".format(characterized))
    print("Uncharacterized hits:        {}".format(uncharacterized))
    print("Unresolved UniProt headers:  {}".format(unresolved))
    print("")
    print("Median length with hit:      {:.1f} aa".format(
        median(hit_lengths)
    ))
    print("Median length without hit:   {:.1f} aa".format(
        median(nohit_lengths)
    ))
    print("")
    print("Pfam among hit proteins:     {} / {}".format(
        hit_pfam,
        n_hits
    ))
    print("Pfam among no-hit proteins:  {} / {}".format(
        nohit_pfam,
        n_nohits
    ))
    print("")
    print("Annotated genes in GFF:      {}".format(annotated_genes))
    print("Annotated mRNAs in GFF:      {}".format(annotated_mrnas))
    print("")
    print("Output GFF:")
    print(args.output_gff)
    print("")
    print("Protein summary:")
    print(args.output_table)
    print("")
    print("========================================")
    print("Completed")
    print("========================================")


if __name__ == "__main__":
    main()
