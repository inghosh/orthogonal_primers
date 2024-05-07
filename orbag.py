from data.params import REF_LIST, PARAMS
from src.primer_interactions import do_primers_interact_full, check_binding_to_ref
from src.primer_generator import generate_random_primers
from itertools import islice
from functools import partial
import multiprocessing as mp


def check_reference_list(query, ref_idx_list):
    statuses = [check_binding_to_ref(query, ref) for ref in ref_idx_list]
    return any(statuses)


def batch(lst, n):
    it = iter(lst)
    return iter(lambda: tuple(islice(it, n)), ())


def grow_barcode_set(
    barcode_list_initial=[], grow_step_size=None, check_references=False, thresh_g=-9
):
    # Generate set of primers
    if grow_step_size is not None:
        primer_set = generate_random_primers(grow_step_size)
    else:
        primer_set = generate_random_primers()

    # Seed barcode_list if empty
    if len(barcode_list_initial) == 0:
        barcode_list = [next(primer_set)]
    else:
        barcode_list = [x for x in barcode_list_initial]

    for count, primer in enumerate(primer_set):
        primer_pass = True
        aligns = False
        for ort_primer_batch in batch(barcode_list, 16):
            if not primer_pass:
                break
            ctx = mp.get_context("fork")
            with ctx.Pool() as pool_obj:
                for result in pool_obj.imap_unordered(
                    partial(
                        do_primers_interact_full,
                        primer1=primer,
                        check_revcomp=PARAMS["check_revcomp"],
                        g_thresh=thresh_g,
                    ),
                    ort_primer_batch,
                ):
                    if result:
                        primer_pass = False

        if primer_pass and check_references:
            aligns = check_reference_list(primer, REF_LIST)
            if not (aligns):
                barcode_list.append(primer)
        elif primer_pass:
            barcode_list.append(primer)
        print(
            f"barcodes in list={len(barcode_list)}",
            f"orthogonal={primer_pass}",
            f"aligns to gDNA={aligns}",
            f"barcodes_tested={count}",
        )

    return barcode_list


def generate_bcs_main():
    # desired_barcodes = PARAMS["desired_barcodes"]
    desired_barcodes = 200
    barcode_bag = []
    while len(barcode_bag) < desired_barcodes:
        print("Generating primers...")
        barcode_bag = grow_barcode_set(barcode_bag)
        with open(f"Output2/{len(barcode_bag)}_barcodes.csv", mode="w") as output_file:
            for barcode in barcode_bag:
                output_file.write(f"{barcode}\n")
        print(len(barcode_bag))
        if len(barcode_bag) > desired_barcodes:
            break
    print(barcode_bag)


if __name__ == "__main__":
    generate_bcs_main()
