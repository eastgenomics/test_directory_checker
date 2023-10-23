import json
from typing import Iterable

import numpy as np
import pandas as pd


# Some Panelapp IDs are not accessible through the Panelapp API either because
# they have been removed or because the panel is in development
UNACCESSIBLE_PANELAPP_IDS = ["481", "1218"]


def parse_td(test_directory, config):
    """Parse rare disease test directory using the config file

    Args:
        test_directory (str): Path to the test directory
        config (dict): Dict containing the data for the config file

    Returns:
        pandas.Dataframe: Dataframe containing the columns of interest
    """

    xls = pd.read_excel(
        test_directory, sheet_name=None, header=config["header_index"]
    )

    for sheet in xls:
        if sheet == config["sheet_of_interest"]:
            data = xls[sheet].loc(axis=1)[
                config["clinical_indication_column_code"],
                config["clinical_indication_column_name"],
                config["panel_column"],
                config["test_method_column"]
            ]

    return data


def parse_hgnc_dump(hgnc_file):
    """ Parse HGNC dump file

    Args:
        hgnc_file (str): Path to the HGNC dump

    Returns:
        pd.Dataframe: Dataframe containing the data in the HGNC dump
    """

    df = pd.read_csv(hgnc_file, delimiter="\t")
    # replace NA by None so that we can handle them
    df_with_none = df.where(pd.notnull(df), None)
    return df_with_none


def parse_genepanels(genepanels):
    """ Parse genepanels file

    Args:
        genepanels (str): Path to the genepanels file

    Returns:
        pd.Dataframe: Dataframe containing the data in the genepanels file
    """

    return pd.read_csv(
        genepanels, delimiter="\t", names=["ci", "panel", "gene"]
    )


def load_config(config):
    """ Load the JSON config for the test directory

    Args:
        config (str): Path to the config file

    Returns:
        dict: Dict containing the data in the JSON config
    """

    config = open(config)
    data = json.load(config)
    config.close()
    return data


def get_all_hgnc_ids_in_target(targets: Iterable, signedoff_panels: dict):
    """ Get the HGNC ids from the panels/genes targets

    Args:
        targets (Iterable): Iterable containing the targets to extract HGNC ids
        from
        signedoff_panels (dict): Dict containing the panelapp ids and the
        corresponding Panelapp panel objects

    Returns:
        set: Set of genes for all the targets
    """

    data = set()

    for target in targets:
        # assume it's a panelapp panel id
        if target.isdigit():
            # some panelapp ids are not accessible through the API because they
            # have been retired or they are in development
            if target in UNACCESSIBLE_PANELAPP_IDS:
                continue

            panel = signedoff_panels[int(target)]
            data.update(
                [gene["hgnc_id"] for gene in panel.get_genes()]
            )

        else:
            # assume it's an HGNC id
            data.add(target)

    return data


def get_genes_from_td_target(
    td_data: pd.DataFrame, signedoff_panels: dict, gene_locus_type: dict
) -> tuple:
    """ Extract the genes from the target columns from the test directory
    either from a Panelapp panel or gene symbols and get their HGNC ids.

    Args:
        td_data (pd.DateFrame): Dataframe containing the test directory data
        signedoff_panels (dict): Dict containing Panelapp IDs as keys and panel
        objects as values
        gene_locus_type (dict): Dict containing the outcome of the gene locus
        type check

    Returns:
        tuple: Tuple containing the set of genes for the given targets and the
        gene locus type dict to get updated
    """

    # Check that the content didn't change
    # get list of HGNC ids for test directory and genepanels
    identified_targets = np.concatenate(
        (
            [
                target
                for sublist in td_data["Identified panels"].to_numpy()
                for target in sublist
            ],
            [
                target
                for sublist in td_data["Identified genes"].to_numpy()
                for target in sublist
            ],
        ), axis=None
    )

    td_genes = set()

    for gene in get_all_hgnc_ids_in_target(
        identified_targets, signedoff_panels
    ):
        if gene in gene_locus_type:
            if gene_locus_type[gene]:
                td_genes.add(gene)

    return td_genes


def filter_out_df(df: pd.DataFrame, **filter_elements) -> pd.DataFrame:
    """ Filter out rows matching values in a dataframe using a dict of data

    Args:
        df (pd.DataFrame): Dataframe to be filtered
        **filter_elements: Kwargs elements to filter with

    Returns:
        pd.DataFrame: Filtered dataframe
    """

    intermediate_data = []

    for key, value in filter_elements.items():
        if value is None:
            df_to_add = df[~df[key].isna()]
        else:
            df_to_add = df[df[key] != value]

        intermediate_data.append(df_to_add)

    filtered_data = pd.concat(intermediate_data).drop_duplicates()

    return filtered_data


def get_locus_status_genes(target_data, signedoff_panels, hgnc_dump):
    """ Extract the genes from the target columns from the test directory
    either from a Panelapp panel or gene symbols and get their HGNC ids.

    Args:
        td_data (pd.DateFrame): Dataframe containing the test directory data
        signedoff_panels (dict): Dict containing Panelapp IDs as keys and panel
        objects as values
        hgnc_dump (pd.DataFrame): Dataframe containing data from the HGNC dump

    Returns:
        tuple: Tuple containing the set of genes for the given targets and the
        gene locus type dict to get updated
    """

    # Check that the content didn't change
    # get list of HGNC ids for test directory and genepanels
    identified_targets = np.concatenate(
        (
            [
                target
                for sublist in target_data["Identified panels"].to_numpy()
                for target in sublist
            ],
            [
                target
                for sublist in target_data["Identified genes"].to_numpy()
                for target in sublist
            ],
        ), axis=None
    )

    gene_locus_type = {}

    for gene in get_all_hgnc_ids_in_target(
        identified_targets, signedoff_panels
    ):
        if gene not in gene_locus_type:
            hgnc_info = hgnc_dump.loc[hgnc_dump["HGNC ID"] == gene]

            # if the gene is TRAC or IGHM, genes that we don't have transcripts
            # for, removed them for the genes for the content comparison
            if gene in ["HGNC:12029", "HGNC:5541"]:
                gene_locus_type[gene] = False
                continue

            # RNA genes and mitochondrial genes are excluded from the
            # genepanels file because we don't have transcripts for
            # them
            if (
                "RNA" in hgnc_info["Locus group"].to_numpy()[0] or
                hgnc_info["Chromosome"].to_numpy()[0] == "mitochondria"
            ):
                gene_locus_type[gene] = False
            else:
                gene_locus_type[gene] = True

    return gene_locus_type
