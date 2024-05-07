import re
import subprocess
import tempfile
from pathlib import Path

from Bio.Seq import Seq
from data.params import RNA_FOLD_PARAMS, PARAMS


def rnacofold_wrapper(primer1, primer2):
    primer_pair_string = f"{primer1}&{primer2}\n"
    parameters_string = f"--paramFile={RNA_FOLD_PARAMS}"
    ps = subprocess.Popen(["echo", primer_pair_string], stdout=subprocess.PIPE)
    folding_output = subprocess.run(
        ["RNAcofold", "--noconv", parameters_string, "--noPS"],
        stdin=ps.stdout,
        capture_output=True,
    )
    ps.wait()
    folding_stdout = folding_output.stdout.decode("utf-8")
    return folding_stdout


def rnaduplex_wrapper(primer, template):
    primer_pair_string = f"{primer}\n{template}\n"
    parameters_string = f"--paramFile={RNA_FOLD_PARAMS}"
    ps = subprocess.Popen(["echo", primer_pair_string], stdout=subprocess.PIPE)
    folding_output = subprocess.run(
        ["RNAduplex", "--noconv", parameters_string],
        stdin=ps.stdout,
        capture_output=True,
    )
    ps.wait()
    folding_stdout = folding_output.stdout.decode("utf-8")
    return folding_stdout


def thresh_g_default(thresh):
    if thresh is None:
        return PARAMS["thresh_g"]
    return thresh


def gibbs_from_fold_output(fold_string):
    gibbs_string = re.split(" ", fold_string)[-1]
    gibbs_value = float(
        "".join([x if x not in " )(\n" else "" for x in list(gibbs_string)])
    )
    return gibbs_value


def does_primer_pair_interact(primer1, primer2, thresh_g=None):
    thresh_g = thresh_g_default(thresh_g)
    folding_result = rnacofold_wrapper(primer1, primer2)
    gibbs_value = gibbs_from_fold_output(folding_result)
    do_primers_interact = gibbs_value < thresh_g
    return do_primers_interact


def does_primer_prime(primer, template, thresh_g=None):
    thresh_g = thresh_g_default(thresh_g)
    folding_result = rnaduplex_wrapper(primer, template)
    gibbs_value = gibbs_from_fold_output(folding_result)
    do_primers_interact = gibbs_value < thresh_g
    return do_primers_interact


def rev_comp_seq(seq):
    fw_seq = Seq(seq)
    rv_seq = fw_seq.reverse_complement()
    return str(rv_seq)


def do_primers_interact_full(primer2, primer1, check_revcomp=True, g_thresh=None):
    if g_thresh is None:
        g_thresh = PARAMS["thresh_g"]
    if does_primer_pair_interact(primer1, primer2, thresh_g=g_thresh):
        return True
    if does_primer_pair_interact(primer2, primer1, thresh_g=g_thresh):
        return True
    if check_revcomp:
        p1_rc = rev_comp_seq(primer1)
        if does_primer_pair_interact(p1_rc, primer2, thresh_g=g_thresh):
            return True
        if does_primer_pair_interact(primer2, p1_rc, thresh_g=g_thresh):
            return True
    return False


def check_binding_to_ref(query_seq, ref_idx_loc):
    with tempfile.TemporaryDirectory() as temp_dir_obj:
        temp_dir = Path(temp_dir_obj)
        read_file = temp_dir / "read_file"
        read_file.write_text(f">_\n{query_seq}")
        output_sam = temp_dir / "output"
        subprocess.run(
            [
                "bowtie2",
                "-f",
                "-k",
                "1",
                "--no-hd",
                "-x",
                ref_idx_loc,
                "-U",
                read_file,
                "-S",
                output_sam,
            ],
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
        )
        sam_line = output_sam.read_text()
        flag = sam_line.split("\t")[1]
    aligned_status = flag != str(4)
    return aligned_status
