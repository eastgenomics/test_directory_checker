import numpy as np
import pandas as pd

from test_directory_checker import utils, identify


def check_target(row: pd.Series, hgnc_dump: pd.DataFrame) -> pd.Series:
    """ Check target and return the identified panels and genes

    Args:
        row (pd.Series): Pandas Series from the test directory
        hgnc_dump (pd.DataFrame): Pandas Dataframe containing the data from the
        HGNC dump file

    Returns:
        pd.Series: Pandas Series containing the new identified panels and
        genes columns
    """

    # stupid weird dash that needs replacing
    target = row["Target/Genes"].replace("–", "-")
    row["Identified panels"], row["Identified genes"] = identify.identify_target(
        target, hgnc_dump
    )
    return row


def check_test_method(row: pd.Series, config: dict) -> pd.Series:
    """ Check the test method from the test directory

    Args:
        row (pd.Series): Pandas Series from the test directory
        config (dict): Dict containing the data from the config file

    Returns:
        pd.Series: Pandas Series containing the potential new test methods and
        removed test methods
    """

    # check for new test methods
    # check for typos
    test_methods_config = config["ngs_test_methods"]
    test_method = [row["Test Method"]]
    diff_potential_new_tm = set(test_method) - set(test_methods_config)
    diff_potential_removed_tm = set(test_methods_config) - set(test_method)
    row["Potential new test methods"] = diff_potential_new_tm
    row["Potential removed test methods"] = diff_potential_removed_tm
    return row


def compare_gp_td(td_data, genepanels_data, hgnc_dump, signedoff_panels):
    identical_ci = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_target", "td_version",
            "td_genes", "removed", "added"
        ]
    )

    replaced_ci = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_target", "td_genes",
            "removed", "added"
        ]
    )

    gene_locus_type = {}

    for gemini_name in genepanels_data["ci"].unique():
        data = {
            "gemini_name": None, "panel": None, "genes": None,
            "td_target": None, "td_genes": None, "removed": None, "added": None
        }
        data["gemini_name"] = gemini_name

        gemini_name_splitted = gemini_name.split("_")
        r_code = gemini_name_splitted[0]

        data_for_r_code = genepanels_data[
            genepanels_data["ci"] == gemini_name
        ]
        genepanels_genes = set(data_for_r_code["gene"].unique())

        data["panel"] = ", ".join(data_for_r_code["panel"].unique())
        data["genes"] = ", ".join(sorted(list(genepanels_genes)))

        if gemini_name.startswith("C"):
            print("'C' clinical indications are bespoke, skipping")
            continue

        td_for_test_id = td_data[
            td_data["Test ID"] == r_code
        ]

        if td_for_test_id.shape[0] == 1:
            # Check that the content didn't change
            # get list of HGNC ids for test directory and genepanels
            identified_targets = np.concatenate(
                (
                    td_for_test_id["Identified panels"].to_numpy()[0],
                    td_for_test_id["Identified genes"].to_numpy()[0]
                ), axis=None
            )

            data["td_target"] = ", ".join(
                td_for_test_id["Target/Genes"].to_numpy()
            )

            td_genes = set()

            for gene in utils.get_all_hgnc_ids_in_target(
                identified_targets, signedoff_panels
            ):
                if gene not in gene_locus_type:
                    hgnc_info = utils.hgnc_query(gene, hgnc_dump)

                    if (
                        "RNA" in hgnc_info["Locus group"].to_numpy()[0] or
                        hgnc_info["Chromosome"].to_numpy()[0] == "mitochondria"
                    ):
                        gene_locus_type[gene] = False
                    else:
                        gene_locus_type[gene] = True
                        td_genes.add(gene)

                else:
                    if gene_locus_type[gene]:
                        td_genes.add(gene)

            data["td_version"] = ", ".join([
                signedoff_panels[int(target)].get_version()
                for target in td_for_test_id["Identified panels"].to_numpy()[0]
                if target != "481"
            ])

            data["td_genes"] = ", ".join(sorted(list(td_genes)))
            removed_genes = genepanels_genes - td_genes
            new_genes = td_genes - genepanels_genes

            if removed_genes:
                data["removed"] = ", ".join(sorted(list(removed_genes)))

            if new_genes:
                data["added"] = ", ".join(sorted(list(new_genes)))

        else:
            print("Check if the test hasn't been replaced by another test")

            td_for_r_code = td_data[
                td_data["Test ID"] == r_code.split(".")[0]
            ]

            if td_for_r_code.shape[0] == 0:
                print("Clinical indication has been removed")
                continue
            elif td_for_r_code.shape[0] == 1:
                print((
                    "Check if clinical indication has been replaced by new "
                    "unique one"
                ))
            elif td_for_r_code.shape[0] >= 2:
                print((
                    "Check if clinical indication has been replaced by one of "
                    "the new ones"
                ))

        identical_ci = identical_ci.append(data, ignore_index=True)

    identical_ci.to_html("ci_existing_in_both.html")
