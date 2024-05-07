import random
from Levenshtein import distance
from orbag import generate_random_dna
from orbag import rev_comp_seq


def homopoly_in_seq(seq):
    if (
        "AAAA" not in seq
        and "CCCC" not in seq
        and "GGGG" not in seq
        and "TTTT" not in seq
    ):
        return False
    else:
        return True


def disllowed_base(seq, disallowed, location):
    if seq[location] == disallowed:
        return True
    else:
        return False


def grow_barcode_set(barcode_list_initial=[], desired_barcodes=96):
    barcode_list = barcode_list_initial
    counter = 0
    while len(barcode_list) < 1:
        counter += 1
        new_seq = generate_random_dna(10)
        if not homopoly_in_seq(new_seq):
            if not disllowed_base(new_seq, "T", 0):
                if not disllowed_base(new_seq, "T", 9):
                    barcode_list.append(new_seq)
                    print(len(barcode_list), new_seq, counter)

    while len(barcode_list) < desired_barcodes:
        counter += 1
        new_seq = generate_random_dna(10)
        if not homopoly_in_seq(new_seq):
            if not disllowed_base(new_seq, "T", 0):
                if not disllowed_base(new_seq, "T", 9):
                    seq_rc = rev_comp_seq(new_seq)
                    match = False
                    for already_here in barcode_list:
                        if (
                            distance(already_here, seq_rc) < 4
                            or distance(already_here, new_seq) < 5
                        ):
                            match = True
                            break
                    if not match:
                        barcode_list.append(new_seq)
                        print(len(barcode_list), new_seq, counter)

    print(barcode_list)

    with open("barcodes.txt", mode="w") as outputfile:
        for barcode in barcode_list:
            outputfile.write(barcode)
            outputfile.write("\n")


if __name__ == "__main__":
    grow_barcode_set()
