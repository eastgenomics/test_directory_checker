import json

import numpy as np
import pandas as pd
import regex


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
