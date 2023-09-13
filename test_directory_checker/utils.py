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
        genepanels, delimiter="\t", columns=["ci", "panel", "gene"]
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


def find_hgnc_id(gene_symbol, hgnc_dump):
    """ Find hgnc id using the hgnc dump

    Args:
        gene_symbol (str): Gene symbol
        hgnc_dump (pd.Dataframe): Hgnc dump dataframe

    Raises:
        Exception: if a panel has escaped previous checks

    Returns:
        str: Hgnc id
    """

    df_res = pd.Series(
        [gene_symbol, None, None, None], index=[
            "Gene symbol", "Selected HGNC ID", "Previous", "Alias"
        ]
    )

    hgnc_id = None

    # pattern is the gene symbol and only the gene symbol
    pattern = fr"^{gene_symbol}$"
    # try and match the gene symbol in the Approved symbol column
    data = hgnc_dump[
        hgnc_dump["Approved symbol"].str.match(pattern, na=False)
    ]["HGNC ID"]
    # prepare dataframes for aliases and previous symbols
    previous_symbols = hgnc_dump[["Previous symbols", "HGNC ID"]]
    alias_symbols = hgnc_dump[["Alias symbols", "HGNC ID"]]

    # match failed, need to use previous and alias symbols
    if len(data.index) == 0:
        if regex.match(r"[A-Z]+[A-Z0-9]+", gene_symbol):
            # previous and alias symbols in hgnc are separated by commas which
            # breaks the pattern matching so i split the appropriate columns
            # creating a df with one column per splitted element
            splitted_previous_symbols = previous_symbols["Previous symbols"].str.split(",", expand=True)
            splitted_alias_symbols = alias_symbols["Alias symbols"].str.split(",", expand=True)
            # copy pasted thing to do the pattern matching for every element of
            # the splitted column
            previous_symbols_match = np.column_stack(
                [
                    splitted_previous_symbols[col].str.strip().str.match(pattern, na=False)
                    for col in splitted_previous_symbols
                ]
            )
            alias_symbols_match = np.column_stack(
                [
                    splitted_alias_symbols[col].str.strip().str.match(pattern, na=False)
                    for col in splitted_alias_symbols
                ]
            )
            # go back to a dataframe using the numpy array to get the matching
            # rows
            df_previous_symbols = previous_symbols.loc[previous_symbols_match.any(axis=1)]
            df_alias_symbols = alias_symbols.loc[alias_symbols_match.any(axis=1)]

            # check the size of the dataframes and do the appropriate action
            if len(df_previous_symbols) == 0 and len(df_alias_symbols) == 0:
                # couldn't find a previous or alias symbol
                return df_res

            elif len(df_previous_symbols) == 1 and len(df_alias_symbols) == 0:
                # found only a previous symbol, return the HGNC id
                hgnc_id = df_previous_symbols["HGNC ID"].to_list()[0]
                df_res.at["Previous"] = True
                df_res.at["Alias"] = False

            elif len(df_previous_symbols) == 0 and len(df_alias_symbols) == 1:
                # found only a alias symbol, return the HGNC id
                hgnc_id = df_alias_symbols["HGNC ID"].to_list()[0]
                df_res.at["Previous"] = False
                df_res.at["Alias"] = True

            elif len(df_previous_symbols) >= 1 and len(df_alias_symbols) >= 1:
                # found previous and alias symbols, cry
                df_res.at["Previous"] = True
                df_res.at["Alias"] = True

            elif len(df_previous_symbols) >= 1:
                df_res.at["Previous"] = True
                df_res.at["Alias"] = False

            elif len(df_alias_symbols) >= 1:
                df_res.at["Previous"] = False
                df_res.at["Alias"] = True

    else:
        hgnc_id = data.iloc[0]

    df_res.at["Selected HGNC ID"] = hgnc_id

    return df_res


def handle_list_panels(panels, hgnc_dump):
    """ Given a list of "panels", get the hgnc ids/rescue comma panelapp panels

    Args:
        panels (list): List of panels
        hgnc_dump (pd.Dataframe): Hgnc dump dataframe

    Returns:
        list: List of hgnc ids/panel
    """

    hgnc_ids = []

    # check if the list only contains genes
    # check that the first element has a gene structure i.e. it's
    # not something weird
    if regex.match(r"[A-Z]+[A-Z0-9]+", panels[0]):
        # first element is a gene, so assume that we're dealing with genes
        for panel in panels:
            # for every element in the list, double check that it is a gene
            # because of test directory unpredictableness
            if regex.match(r"[A-Z]+[A-Z0-9]+", panel):
                hgnc_id = find_hgnc_id(panel, hgnc_dump)

                if hgnc_id:
                    hgnc_ids.append(hgnc_id)
                else:
                    # there are instances where the targets can include "and"
                    # at the end of a gene list.
                    # try and catch those instances
                    if "and" in panel:
                        attempt_rescue_gene = [
                            gene.strip() for gene in panel.split("and")
                        ]

                        # check if there were 2 genes from each side of the 'and'
                        if len(attempt_rescue_gene) >= 2:
                            for gene in attempt_rescue_gene:
                                rescued_gene = find_hgnc_id(gene, hgnc_dump)

                                if rescued_gene:
                                    hgnc_ids.append(rescued_gene)
                                else:
                                    hgnc_ids.append(None)
                    else:
                        # we didn't manage to find a HGNC id for the gene
                        # symbol
                        print(
                            f"Couldn't find a HGNC id for '{panel}', please "
                            "check manually"
                        )
                        hgnc_ids.append(None)
            else:
                # that element of the list is not a gene
                print(
                    f"This element '{panel}' was not detected as being a gene."
                    " Please check manually"
                )
                hgnc_ids.append(None)

        return hgnc_ids
    else:
        # do regex to see if every element in the list is a panel
        matches = [
            regex.match(r"[A-Za-z0-9-()\ ,]*\([0-9&\ ]+\)", panel)
            for panel in panels
        ]

        if all(matches):
            # all the elements are panelapp panels
            return extract_panelapp_id(panels)
        else:
            # working on it, i realised that trying to rescue the potential
            # panelapp panels is not trivial at all i.e. what if the panel has
            # multiple commas
            # so i think it's safer to manual rescue them
            print(
                f"Potential panelapp panels with commas '{panels}'"
            )
            return None


def extract_panelapp_id(panel_match):
    """ Extract the panelapp id from the target in the test directory

    Args:
        panels (iter): Iterable containing the panels to look at

    Returns:
        list: List of panelapp ids
    """

    panelapp_id_match = regex.match(
        r"(?P<panelapp_id>\([0-9]+\))", panel_match
    )

    if panelapp_id_match:
        cleaned_panelapp_id = panelapp_id_match.group(
            "panelapp_id"
        ).replace("(", "").replace(")", "")

        return cleaned_panelapp_id

    return None


def write_json(dict_data):
    with open("targets.json", "w") as f:
        json.dump(dict_data, f, indent=2)


def write_test_methods(new_test_methods, removed_test_methods):
    with open("potential_new_test_methods.txt", "w") as f:
        for tm in sorted(new_test_methods):
            f.write(f"{tm}\n")

    with open("potential_removed_test_methods.txt", "w") as f:
        for tm in sorted(removed_test_methods):
            f.write(f"{tm}\n")
