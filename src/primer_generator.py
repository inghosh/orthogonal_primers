import random
from pathlib import Path
import re
import subprocess
import tempfile

from data.params import PARAMS


def generate_random_dna(length=1000):
    """Generates a stretch of random DNA"""

    random_dna = "".join(random.choice("CGTA") for _ in range(length))
    return random_dna


BLANK_SETTINGS_FILE = Path("../data/p3_settings_blank.txt")


def parse_boulder_IO(input):
    """Converts a boulder IO formatted string into a dict of entries"""

    output = {}
    for line in re.split("\n", input):
        entry = re.split("=", line)
        if entry[0] != "":
            output[entry[0]] = entry[-1]
    return output


def primer_list_from_primer3(primer3_output):
    """parses the Boulder IO output from primer3 to generate a list of primers"""

    parsed_data = parse_boulder_IO(primer3_output)
    for key, value in parsed_data.items():
        if re.search("_SEQUENCE", key) is not None:
            yield value


def generate_primers_from_dna(
    input_seq, primer_count_requested=None, tm_opt=None, end_stability=None
):
    """Given a stretch of DNA, yields several conforming primer sequences"""

    if primer_count_requested is None:
        primer_count_requested = int(len(input_seq) / PARAMS["rand_seq_per_prim"])

    if tm_opt is None:
        tm_opt = PARAMS["opt_tm"]

    tm_min = tm_opt - 3
    tm_max = tm_opt + 3

    if end_stability is None:
        end_stability = PARAMS["end_stability"]

    template_line = f"SEQUENCE_TEMPLATE={input_seq}\n"
    primer_num_line = f"PRIMER_NUM_RETURN={primer_count_requested}\n"

    tm_opt_line = f"PRIMER_OPT_TM={tm_opt}\n"
    tm_min_line = f"PRIMER_MIN_TM={tm_min}\n"
    tm_max_line = f"PRIMER_MAX_TM={tm_max}\n"

    end_stability_line = f"PRIMER_MAX_END_STABILITY={end_stability}\n"

    with tempfile.NamedTemporaryFile(mode="w+t") as temp_file_obj:
        temp_file_loc = Path(temp_file_obj.name)
        temp_file_obj.write(template_line)
        temp_file_obj.write(primer_num_line)
        temp_file_obj.write(tm_opt_line + tm_min_line + tm_max_line)
        temp_file_obj.write(end_stability_line)
        temp_file_obj.write(BLANK_SETTINGS_FILE.read_text())
        temp_file_obj.seek(0)
        primer3_output = subprocess.run(
            ["primer3_core", temp_file_loc.resolve()], capture_output=True
        )
        primer3_stdout = primer3_output.stdout.decode("utf-8")
        primer_list = primer_list_from_primer3(primer3_stdout)

    for primer in primer_list:
        yield primer


def generate_random_primers(n=None, rand_seq_len=None):
    """generates a bunch of primers, making 'n' at a time"""
    if n is None:
        n = PARAMS["primers_gen_per_set"]
    if rand_seq_len is None:
        rand_seq_len = n * PARAMS["rand_seq_per_prim"]

    primer_gen = generate_primers_from_dna(generate_random_dna(rand_seq_len), n)
    for primer in primer_gen:
        yield primer
