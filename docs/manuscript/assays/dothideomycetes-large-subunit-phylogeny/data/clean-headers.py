import re
from pathlib import Path

from Bio import SeqIO


def __main(
    input_file: Path,
    output_file: Path,
) -> None:
    if output_file.exists():
        output_file.unlink()

    with input_file.open("r") as f:
        for record in SeqIO.parse(f, "fasta"):
            record.id = record.id.split(".")[0]
            record.description = re.sub(
                r"[^a-zA-Z0-9_]",
                "",
                " ".join(record.description.split(" ")[1:5]).replace(" ", "_"),
            )

            with output_file.open("a+") as f:
                f.write(f">{record.id}_{record.description}\n{record.seq}\n")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Clean headers of FASTA file downloaded from GenBank"
    )

    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    __main(input_file=Path(args.input), output_file=Path(args.output))
