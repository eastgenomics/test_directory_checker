import regex

import numpy as np
import pandas as pd


def identify_target(target: str, hgnc_dump: pd.DataFrame) -> list:
    """ Identify the target as gene or panel using regex

    Args:
        target (str): String for the target extracted from the test directory
        hgnc_dump (pd.DataFrame): Pandas Dataframe containing the HGNC data

    Returns:
        list: List of 2 elements containing the identified panels and the
        identified genes
    """

    panels = []
    genes = []

    potential_panel_targets = regex.findall(r"\([0-9&\ ]+\)", target)
    potential_gene_targets = regex.findall(r"[A-Z]+[A-Z0-9\-]+", target)

    # regex to identify panelapp panels
    if potential_panel_targets:
        for potential_panel in potential_panel_targets:
            cleaned_panelapp_id = potential_panel.replace(
                "(", "").replace(")", "")
            panels.append(cleaned_panelapp_id)

    # regex to identify gene symbol
    elif potential_gene_targets:
        for potential_gene in potential_gene_targets:
            hgnc_id_data = find_hgnc_id(potential_gene, hgnc_dump)

            if hgnc_id_data["HGNC ID"]:
                genes.append(hgnc_id_data["HGNC ID"])

    return panels, genes


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
            "Gene symbol", "HGNC ID", "Previous", "Alias"
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

    df_res.at["HGNC ID"] = hgnc_id

    return df_res
