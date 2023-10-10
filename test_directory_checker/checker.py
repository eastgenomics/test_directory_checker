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
    row["Potential new test methods"] = diff_potential_new_tm
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

    identical_tests = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target", "td_version",
            "td_genes", "removed", "added"
        ]
    )

    removed_tests = pd.DataFrame([], columns=["gemini_name", "panel", "genes"])

    replaced_tests = pd.DataFrame(
        [],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target", "td_genes",
            "removed", "added"
        ]
    )

    # this dict will contain genes and whether this gene is captured by our
    # analysis downstream due to its locus type i.e. RNA or mitochondrial
    # it is used to automatically skip those RNA and mitochondrial genes.
    gene_locus_type = {}

    # go through every test ID in the genepanels file
    for gemini_name in genepanels_data["ci"].unique():
        data = {
            "gemini_name": None, "panel": None, "genes": None, "td_ci": None,
            "td_target": None, "td_genes": None, "removed": None, "added": None
        }
        data["gemini_name"] = gemini_name

        gemini_name_splitted = gemini_name.split("_")
        r_code = gemini_name_splitted[0]

        # get subset using the gemini name
        data_for_r_code = genepanels_data[
            genepanels_data["ci"] == gemini_name
        ]
        # get the genes for that test ID
        genepanels_genes = set(data_for_r_code["gene"].unique())

        data["panel"] = ", ".join(data_for_r_code["panel"].unique())
        data["genes"] = ", ".join(sorted(list(genepanels_genes)))

        if gemini_name.startswith("C"):
            print("'C' clinical indications are bespoke, skipping")
            continue

        # filter td data using the r-code
        td_for_test_id = td_data[
            td_data["Test ID"] == r_code
        ]

        if td_for_test_id.shape[0] == 1:
            # found test id in test directory
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

            identical_tests = identical_tests.append(data, ignore_index=True)

        else:
            # didn't find the test ID, use clinical indication ID to find
            # equivalence
            td_for_r_code = td_data[
                td_data["Test ID"].str.contains(r_code.split(".")[0])
            ]

            if td_for_r_code.shape[0] == 0:
                # clinical indication has been removed
                removed_tests = removed_tests.append(data, ignore_index=True)
                removed_tests = removed_tests.reindex(
                    columns=["gemini_name", "panel", "genes"]
                )

            elif td_for_r_code.shape[0] == 1:
                # check if that new test code replaces the old one by looking
                # at the gene content
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

                replaced_tests = replaced_tests.append(data, ignore_index=True)

            elif td_for_r_code.shape[0] >= 2:
                # loop through those tests and check if one of them replaces
                # the old one
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

                    replaced_tests = replaced_tests.append(
                        data, ignore_index=True
                    )

    identical_tests = identical_tests.reindex(
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    replaced_tests = replaced_tests.reindex(
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    return identical_tests, removed_tests, replaced_tests


def find_new_clinical_indications(
    td_data: pd.DataFrame, genepanels_df: pd.DataFrame
) -> pd.DataFrame:
    """ Find new clinical indications i.e. present in the test directory but
    not in the genepanels file

    Args:
        td_data (pd.DataFrame): Dataframe with test directory data
        genepanels_df (pd.DataFrame): Dataframe with genepanels data

    Returns:
        pd.DataFrame: Dataframe containing the new tests from the test
        directory
    """

    # extract the r codes from the genepanels file
    genepanels_rcodes = [ci.split("_")[0] for ci in genepanels_df["ci"].values]
    return td_data[~td_data["Test ID"].isin(genepanels_rcodes)]
