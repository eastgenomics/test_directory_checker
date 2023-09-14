import regex

from test_directory_checker import utils


def identify_target(target, hgnc_dump):

    res = {}
    res.setdefault("panels", [])
    res.setdefault("genes", [])

    potential_panel_targets = regex.findall(r"\([0-9&\ ]+\)", target)
    potential_gene_targets = regex.findall(r"[A-Z]+[A-Z0-9]+", target)

    # regex to identify panelapp panels
    if potential_panel_targets:
        for potential_panel in potential_panel_targets:
            cleaned_panelapp_id = potential_panel.replace(
                "(", "").replace(")", "")
            res["panels"].append(cleaned_panelapp_id)

    # regex to identify gene symbol
    if potential_gene_targets:
        for potential_gene in potential_gene_targets:
            hgnc_id_data = utils.find_hgnc_id(potential_gene, hgnc_dump)

            if hgnc_id_data["HGNC ID"]:
                res["genes"].append(hgnc_id_data["HGNC ID"])

    return res
