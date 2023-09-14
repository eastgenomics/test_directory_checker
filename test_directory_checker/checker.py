import regex

from test_directory_checker import utils, identify


def check_target(row, hgnc_dump):
    # stupid weird dash that needs replacing
    target = row["Target/Genes"].replace("–", "-")
    targets = identify.identify_target(target, hgnc_dump)
    row["Identified targets"] = targets
    return row


def check_test_method(row, config):
    # check for new test methods
    # check for typos
    test_methods_config = config["ngs_test_methods"]
    test_method = [row["Test Method"]]
    diff_potential_new_tm = set(test_method) - set(test_methods_config)
    diff_potential_removed_tm = set(test_methods_config) - set(test_method)
    row["Potential new test methods"] = diff_potential_new_tm
    row["Potential removed test methods"] = diff_potential_removed_tm
    return row


def compare_gp_td(td_data, genepanels_data, hgnc_dump):
    for gemini_name in genepanels_data["ci"].unique():
        r_code, ci_name = gemini_name.split("_")
        data_for_r_code = genepanels_data[
            genepanels_data["ci"] == gemini_name
        ]

        td_for_test_id = td_data[
            td_data["Test ID"] == r_code
        ]

        if td_for_test_id.shape[0] == 1:
            print("Check that the content didn't change")
            targets = td_for_test_id["Target/Genes"]
            check_targets(targets, hgnc_dump)
        else:
            print("Check if the test hasn't been replaced by another test")

            td_for_r_code = td_data[
                td_data["Test ID"] == r_code.split(".")[0]
            ]

            if td_for_r_code.shape[0] == 0:
                print("Clinical indication has been removed")
                continue
            elif td_for_r_code.shape[0] == 1:
                print((
                    "Check if clinical indication has been replaced by new "
                    "unique one"
                ))
            elif td_for_r_code.shape[0] >= 2:
                print((
                    "Check if clinical indication has been replaced by one of "
                    "the new ones"
                ))

        for panel in data_for_r_code["panels"].unique():
            genes = data_for_r_code[
                data_for_r_code["panels"] == panel
            ]["genes"].unique()
