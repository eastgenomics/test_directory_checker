import regex

from test_directory_checker import utils


def check_targets(targets, hgnc_dump):
    panels_list = {}

    # check that the panel ID matches the panel name
    for i, target in targets.items():
        # stupid weird dash that needs replacing
        panels = target.replace("–", "-")
        panels_comma = [p.strip() for p in panels.split(",")]
        panels_semicolon = [p.strip() for p in panels.split(";")]

        if panels_comma == panels_semicolon:
            # regex to identify panelapp panels
            if regex.match(
                r"[A-Za-z0-9-()\ ]*\([0-9&\ ]+\)", panels_comma[0]
            ):
                res = utils.extract_panelapp_id(panels_comma)
                panels_list.setdefault("panelapp", {})
                panels_list["panelapp"].setdefault(list(res.keys())[0], []).extend(
                    res[list(res.keys())[0]]
                )
                continue

            # regex to identify gene symbol
            if regex.match(r"[A-Z]+[A-Z0-9]+", panels_comma[0]):
                hgnc_id_data = utils.find_hgnc_id(panels_comma[0], hgnc_dump)
                panels_list.setdefault("gene", {})
                panels_list["gene"].setdefault(list(hgnc_id_data.keys())[0], []).append(
                    hgnc_id_data[list(hgnc_id_data.keys())[0]]
                )
                continue

            # regex to identify the rest
            if regex.match(r"[A-Za-z\ ]", panels_comma[0]):
                panels_list.setdefault("other", []).append(panels_comma)
                continue

        else:
            if len(panels_comma) == 1:
                # try and rescue some panelapp panels
                if regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_comma[0]
                ):
                    res = utils.extract_panelapp_id(panels_comma)
                    panels_list.setdefault("panelapp", {})
                    panels_list["panelapp"].setdefault(list(res.keys())[0], []).extend(
                        res[list(res.keys())[0]]
                    )
                    continue

                else:
                    # assume that we have lists of genes using semicolon
                    pass
                    # print("assume lists of gene semicolon", panels_comma)

            elif len(panels_comma) >= 2:
                cleaned_panels = utils.handle_list_panels(
                    panels_comma, hgnc_dump
                )

                if cleaned_panels:
                    for ele in cleaned_panels:
                        panels_list.setdefault("gene", {})
                        panels_list["gene"].setdefault(
                            list(ele.keys())[0], []
                        ).append(
                            ele[list(ele.keys())[0]]
                        )
                    continue

            if len(panels_semicolon) == 1:
                # try and rescue some panelapp panels
                if regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_semicolon[0]
                ):
                    res = utils.extract_panelapp_id(panels_semicolon)
                    panels_list.setdefault("panelapp", {})
                    panels_list["panelapp"].setdefault(list(res.keys())[0], []).extend(
                        res[list(res.keys())[0]]
                    )
                    continue

                else:
                    # assume that we have lists of genes not using comma
                    # print("assume lists of gene comma", panels_semicolon)
                    pass

            elif len(panels_semicolon) >= 2:
                cleaned_panels = utils.handle_list_panels(
                    panels_semicolon, hgnc_dump
                )

                if cleaned_panels:
                    for ele in cleaned_panels:
                        panels_list.setdefault("gene", {})
                        panels_list["gene"].setdefault(
                            list(ele.keys())[0], []
                        ).append(
                            ele[list(ele.keys())[0]]
                        )
                    continue

    return panels_list


def check_test_methods(test_methods, config):
    # check for new test methods
    # check for typos
    test_methods_td = sorted(list(set(test_methods.values())))
    test_methods_config = config["ngs_test_methods"]
    diff_potential_new_tm = set(test_methods_td) - set(test_methods_config)
    diff_potential_removed_tm = set(test_methods_config) - set(test_methods_td)
    return diff_potential_new_tm, diff_potential_removed_tm
