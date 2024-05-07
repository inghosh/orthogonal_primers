import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def check_bowtie2():
    try:
        subprocess.check_output(["which", "bowtie2"])
    except OSError:
        raise RuntimeError("bowtie2 not found")


def bowtie2_align_generator(query: SeqRecord, subject: SeqRecord):
    """aligns a fasta subject sequence to a query sequence using bowtie2 and generates sam lines"""
    check_bowtie2()
    with tempfile.TemporaryDirectory() as temp_dir_obj:
        temp_dir = Path(temp_dir_obj)
        fasta_reference = temp_dir / "query.fasta"
        SeqIO.write(query, fasta_reference, "fasta")
        index_reference = temp_dir / "index"
        subprocess.run(
            ["bowtie2-build", fasta_reference, index_reference],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        subject_sequence = temp_dir / "subject"
        SeqIO.write(subject, subject_sequence, "fasta")
        output_sam = temp_dir / "output"
        subprocess.run(
            [
                "bowtie2",
                "-f",
                "-a",
                "--no-hd",
                "--very-sensitive-local",
                "-x",
                index_reference,
                "-U",
                subject_sequence,
                "-S",
                output_sam,
            ],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        with open(output_sam) as output:
            lines = output.readlines()
            for line in lines:
                yield line
