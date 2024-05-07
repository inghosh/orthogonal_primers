from pathlib import Path

RNA_FOLD_PARAMS = Path("dna_mathews2004.par")
ECOLI_GDNA_IDX = Path("data/gDNA/ecoli/ec_gdna_idx")
REF_LIST = [ECOLI_GDNA_IDX]
PARAMS = {
    "desired_barcodes": 500,
    "opt_tm": 55,
    "thresh_g": -10,
    "check_revcomp": True,
    "end_stability": 8,
    "check_gdna": False,
    "primers_gen_per_set": 50,
    "rand_seq_per_prim": 100,
}
