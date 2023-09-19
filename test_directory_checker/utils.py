import json
from typing import Iterable

import pandas as pd


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


def hgnc_query(hgnc_id, hgnc_dump):
    row = hgnc_dump.loc[hgnc_dump["HGNC ID"] == hgnc_id]
    return row


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
        if target.isdigit():
            # assume it's a panelapp panel id
            # 481 has been merged with 480
            if target == "481":
                continue

            panel = signedoff_panels[int(target)]
            data.update(
                [gene["hgnc_id"] for gene in panel.get_genes()]
            )

        else:
            # assume it's an HGNC id
            data.add(target)

    return data
