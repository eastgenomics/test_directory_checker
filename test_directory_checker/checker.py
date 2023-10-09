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


def compare_gp_td(
    td_data: pd.DataFrame, genepanels_data: pd.DataFrame,
    hgnc_dump: pd.DataFrame, signedoff_panels: dict 
) -> tuple:
    """ Compare the test directory data and the genepanels data.
    The code will look for test IDs and will compare the content resulting in 3
    test cases:
    - Found test ID in both test directory and genepanels --> compare content
    - Didn't find test ID in both, look for clinical indication ID
        - Didn't find clinical indication ID --> test has been removed
        - Found clinical indication ID --> compare content of test IDs to check
          if we can reproduce the same target content in genepanels

    Args:
        td_data (pd.DataFrame): Dataframe with the test directory data
        genepanels_data (pd.DataFrame): Dataframe with the genepanels data
        hgnc_dump (pd.DataFrame): Dataframe with the HGNC data
        signedoff_panels (dict): Dict containing the signedoff panels from
        Panelapp with the key being the Panelapp ID and the value being a Panel
        object

    Returns:
        tuple: Tuple of 2 elements containing the identical/replaced test IDs
        and the results of the content comparison.
    """

    identical_ci = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target", "td_version",
            "td_genes", "removed", "added"
        ]
    )

    replaced_ci = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target", "td_genes",
            "removed", "added"
        ]
    )

    gene_locus_type = {}

    for gemini_name in genepanels_data["ci"].unique():
        data = {
            "gemini_name": None, "panel": None, "genes": None, "td_ci": None,
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
            td_genes, gene_locus_type_update = utils.get_genes_from_td_target(
                td_for_test_id, signedoff_panels, hgnc_dump, gene_locus_type
            )

            gene_locus_type = {**gene_locus_type, **gene_locus_type_update}

            data["td_ci"] = ", ".join(
                td_for_test_id["Test ID"].to_numpy()
            )
            data["td_target"] = ", ".join(
                td_for_test_id["Target/Genes"].to_numpy()
            )
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

            identical_ci = identical_ci.append(data, ignore_index=True)

        else:
            print("Check if the test hasn't been replaced by another test")

            td_for_r_code = td_data[
                td_data["Test ID"].str.contains(r_code.split(".")[0])
            ]

            if td_for_r_code.shape[0] == 0:
                print("Clinical indication has been removed")
                continue
            elif td_for_r_code.shape[0] == 1:
                print((
                    "Check if clinical indication has been replaced by new "
                    "unique one"
                ))
                td_genes, gene_locus_type_update = utils.get_genes_from_td_target(
                    td_for_r_code, signedoff_panels, hgnc_dump, gene_locus_type
                )

                gene_locus_type = {**gene_locus_type, **gene_locus_type_update}

                data["td_ci"] = ", ".join(
                    td_for_r_code["Test ID"].to_numpy()
                )
                data["td_target"] = ", ".join(
                    td_for_r_code["Target/Genes"].to_numpy()
                )
                data["td_version"] = ", ".join([
                    signedoff_panels[int(target)].get_version()
                    for target in td_for_r_code["Identified panels"].to_numpy()[0]
                    if target != "481"
                ])
                data["td_genes"] = ", ".join(sorted(list(td_genes)))

                removed_genes = genepanels_genes - td_genes
                new_genes = td_genes - genepanels_genes

                if removed_genes:
                    data["removed"] = ", ".join(sorted(list(removed_genes)))

                if new_genes:
                    data["added"] = ", ".join(sorted(list(new_genes)))

                replaced_ci = replaced_ci.append(data, ignore_index=True)

            elif td_for_r_code.shape[0] >= 2:
                print((
                    "Check if clinical indication has been replaced by one of "
                    "the new ones"
                ))

                for i, row in td_for_r_code.iterrows():
                    df = row.to_frame().T

                    td_genes, gene_locus_type_update = utils.get_genes_from_td_target(
                        df, signedoff_panels, hgnc_dump, gene_locus_type
                    )
                    gene_locus_type = {
                        **gene_locus_type, **gene_locus_type_update
                    }

                    data["td_ci"] = ", ".join(
                        df["Test ID"].to_numpy()
                    )
                    data["td_target"] = ", ".join(
                        df["Target/Genes"].to_numpy()
                    )
                    data["td_version"] = ", ".join([
                        signedoff_panels[int(target)].get_version()
                        for target in df["Identified panels"].to_numpy()[0]
                        if target != "481"
                    ])
                    data["td_genes"] = ", ".join(sorted(list(td_genes)))

                    removed_genes = genepanels_genes - td_genes
                    new_genes = td_genes - genepanels_genes

                    if removed_genes:
                        data["removed"] = ", ".join(sorted(list(removed_genes)))

                    if new_genes:
                        data["added"] = ", ".join(sorted(list(new_genes)))

                    replaced_ci = replaced_ci.append(data, ignore_index=True)

    identical_ci = identical_ci.reindex(
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    replaced_ci = replaced_ci.reindex(
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    return identical_ci, replaced_ci


def find_new_clinical_indications(td_data, genepanels_df):
    new_cis = pd.DataFrame(
        [],
        columns=[
            "td_ci", "td_target", "identified_panels", "identified_genes"
        ]
    )

    for test_code in td_data["Test ID"].unique():
        data = {
            "td_ci": None, "td_target": None, "identified_panels": None,
            "identified_genes": None
        }
        td_for_test_id = td_data[td_data["Test ID"] == test_code]
        genepanels_data = genepanels_df[
            genepanels_df["ci"].str.contains(test_code)
        ]

        if genepanels_data.shape[0] == 0:
            data["td_ci"] = ", ".join(
                td_for_test_id["Test ID"].to_numpy()
            )
            data["td_target"] = ", ".join(
                td_for_test_id["Target/Genes"].to_numpy()
            )
            data["identified_panels"] = ", ".join(
                [
                    panel
                    for panel in td_for_test_id["Identified panels"].to_numpy()[0]
                ]
            )
            data["identified_genes"] = ", ".join(
                [
                    gene
                    for gene in td_for_test_id["Identified genes"].to_numpy()[0]
                ]
            )
            new_cis = new_cis.append(data, ignore_index=True)

    return new_cis
