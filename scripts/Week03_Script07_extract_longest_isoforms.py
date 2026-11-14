#!/usr/bin/env python3

from pathlib import Path

FINALDIR = Path("output/MAKER/final")

INPUTS = {
    "protein": FINALDIR / "assembly.all.maker.proteins.fasta.renamed.filtered.fasta",
    "transcript": FINALDIR / "assembly.all.maker.transcripts.fasta.renamed.filtered.fasta",
}

OUTPUTS = {
    "protein": FINALDIR / "assembly.all.maker.proteins.fasta.renamed.longest.fasta",
    "transcript": FINALDIR / "assembly.all.maker.transcripts.fasta.renamed.longest.fasta",
}


def read_fasta(path):
    header = None
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)

                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if header is not None:
            yield header, "".join(sequence)


def gene_id_from_header(header):
    seq_id = header.split()[0]

    if "-R" not in seq_id:
        raise ValueError(
            "Could not determine gene ID from sequence ID: {}".format(seq_id)
        )

    return seq_id.split("-R", 1)[0]


def write_longest(input_path, output_path, label):
    longest = {}
    total_records = 0

    for header, sequence in read_fasta(input_path):
        total_records += 1
        gene_id = gene_id_from_header(header)

        if gene_id not in longest:
            longest[gene_id] = (header, sequence)
        elif len(sequence) > len(longest[gene_id][1]):
            longest[gene_id] = (header, sequence)

    with output_path.open("w") as out:
        for gene_id in sorted(longest):
            header, sequence = longest[gene_id]

            out.write(">{}\n".format(header))

            for i in range(0, len(sequence), 60):
                out.write(sequence[i:i + 60] + "\n")

    print("{}:".format(label))
    print("  Input records:    {}".format(total_records))
    print("  Genes retained:   {}".format(len(longest)))
    print("  Output:           {}".format(output_path))
    print()


def main():
    print("========================================")
    print("Week 3 - Extract longest isoforms")
    print("========================================")
    print()

    for key in ("protein", "transcript"):
        if not INPUTS[key].exists():
            raise FileNotFoundError(str(INPUTS[key]))

        write_longest(
            INPUTS[key],
            OUTPUTS[key],
            key.capitalize()
        )

    print("========================================")


if __name__ == "__main__":
    main()
