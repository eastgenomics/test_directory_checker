from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.schema import MetaData

from test_directory_checker import utils, identify


def check_target(
    test_directory_row: pd.Series, hgnc_dump: pd.DataFrame
) -> pd.Series:
    """ Check the target column and return the identified panels and genes

    Args:
        test_directory_row (pd.Series): Pandas Series from the test directory
        hgnc_dump (pd.DataFrame): Pandas Dataframe containing the data from the
        HGNC dump file

    Returns:
        pd.Series: Pandas Series containing the new identified panels and
        genes columns
    """

    (
        test_directory_row["Identified panels"],
        test_directory_row["Identified genes"]
    ) = identify.identify_target(test_directory_row["Target/Genes"], hgnc_dump)

    return test_directory_row


def check_test_method(
    test_directory_row: pd.Series, config: dict
) -> pd.Series:
    """ Check the test method from the test directory by looking at the list of
    covered test methods in the test directory parser config file

    Args:
        test_directory_row (pd.Series): Pandas Series from the test directory
        config (dict): Dict containing the data from the config file

    Returns:
        pd.Series: Pandas Series containing the potential new test methods and
        removed test methods
    """

    # check for new test methods
    # check for typos
    test_methods_config = config["ngs_test_methods"]
    test_method = [test_directory_row["Test Method"]]
    diff_potential_new_tm = set(test_method) - set(test_methods_config)
    test_directory_row["Potential new test methods"] = ", ".join(
        sorted(list(diff_potential_new_tm))
    )
    return test_directory_row


def compare_gp_td(
    td_data: pd.DataFrame, genepanels_data: pd.DataFrame,
    signedoff_panels: dict, gene_locus_type: dict, blacklist_config: dict
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
        signedoff_panels (dict): Dict containing the signedoff panels from
        Panelapp with the key being the Panelapp ID and the value being a Panel
        object
        gene_locus_type (dict): Dict containing the genes and whether we
        capture them according to their locus type
        blacklist_config (dict): Dict containing blacklisted panel IDs or genes

    Returns:
        tuple: Tuple of 2 elements containing the identical/replaced test IDs
        and the results of the content comparison.
    """

    identical_tests_data = []
    removed_tests_data = []
    replaced_tests_data = []

    # go through every test ID in the genepanels file
    for gemini_name in genepanels_data["ci"].unique():
        if gemini_name.startswith("C"):
            print("'C' clinical indications are bespoke, skipping")
            continue

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

        # filter td data using the r-code
        td_for_test_id = td_data[
            td_data["Test ID"] == r_code
        ]

        if td_for_test_id.shape[0] == 1:
            # found test id in test directory
            data_for_test_id = utils.format_td_data(
                td_for_test_id, genepanels_genes, signedoff_panels,
                gene_locus_type, blacklist_config
            )

            data = {**data, **data_for_test_id}
            identical_tests_data.append(data)

        else:
            # didn't find the test ID, use clinical indication ID to find
            # equivalence
            td_for_r_code = td_data[
                td_data["Test ID"].str.contains(r_code.split(".")[0])
            ]

            if td_for_r_code.shape[0] == 0:
                # clinical indication has been removed
                copy_data = data.copy()
                removed_tests_data.append(copy_data)

            elif td_for_r_code.shape[0] == 1:
                # check if that new test code replaces the old one by looking
                # at the gene content
                data_for_r_code = utils.format_td_data(
                    td_for_r_code, genepanels_genes, signedoff_panels,
                    gene_locus_type, blacklist_config
                )

                data = {**data, **data_for_r_code}
                replaced_tests_data.append(data)

            elif td_for_r_code.shape[0] >= 2:
                # loop through those tests and check if one of them replaces
                # the old one
                for i, row in td_for_r_code.iterrows():
                    retain_data = data

                    # transpose the dataframe so that indexes becomes columns
                    df_for_row = row.to_frame().T

                    data_for_row = utils.format_td_data(
                        df_for_row, genepanels_genes, signedoff_panels,
                        gene_locus_type, blacklist_config
                    )

                    retain_data = {**retain_data, **data_for_row}
                    replaced_tests_data.append(retain_data)

    identical_tests_df = pd.DataFrame(
        identical_tests_data,
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    removed_tests_df = pd.DataFrame(
        removed_tests_data, columns=["gemini_name", "panel", "genes"]
    )

    replaced_tests_df = pd.DataFrame(
        replaced_tests_data,
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    return (
        identical_tests_df, removed_tests_df, replaced_tests_df
    )


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


def check_if_genes_present_in_db(
    username: str, pwd: str, db: str, db_type: str, genes: set
):
    """ Check if the genes in a Series are present in a given database as well
    as checking if a gene is in the database if it does have a clinical
    transcript

    Args:
        username (str): Username to the database
        pwd (str): Corresponding password to the username
        db (str): Database name
        db_type (str): Database type, either MySQL or SQLite
        genes (set): Set containing the gene data

    Raises:
        e: _description_

    Returns:
        _type_: _description_
    """

    try:
        if db_type == "mysql":
            db = create_engine(f"mysql://{username}:{pwd}@localhost/{db}")

        elif db_type == "sqlite":
            assert Path(f"{db}").exists(), f"'{db}' doesn't exist"
            db = create_engine(f"sqlite:///{db}")

    except Exception as e:
        raise e
    else:
        meta = MetaData()
        meta.reflect(bind=db)
        Session = sessionmaker(bind=db)
        session = Session()

    gene_tb = meta.tables["gene"]
    g2t_tb = meta.tables["genes2transcripts"]

    data = []

    for gene in genes:
        # query the database using the HGNC id and join the gene and g2t tables
        query = session.query(gene_tb.c.hgnc_id, g2t_tb.c.clinical_transcript)\
            .join(g2t_tb)\
            .filter(gene_tb.c.hgnc_id == gene).all()

        has_clinical_tx = False
        presence_in_db = False

        if query:
            presence_in_db = True

            # for every hgnc id/transcript couple, check if the clinical
            # transcript status == 1
            # i.e. the transcript is the clinical transcript
            for hgnc_id, clinical_transcript in query:
                if clinical_transcript == 1:
                    has_clinical_tx = True

        data.append([gene, presence_in_db, has_clinical_tx])

    data = sorted(data, key=lambda x: x[0])

    df = pd.DataFrame(
        data, columns=["gene", "presence_in_db", "has_clinical_transcript"]
    )

    return df
