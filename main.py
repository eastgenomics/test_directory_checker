import argparse

from panelapp import queries

from test_directory_checker import checker, utils


def main(args):
    config = utils.load_config("configs/column_config.json")
    td_config = utils.load_config(args["config"])
    td_data = utils.parse_td(args["test_directory"], config)
    hgnc_data = utils.parse_hgnc_dump(args["hgnc_dump"])
    genepanels_data = utils.parse_genepanels(args["genepanels"])
    signedoff_panels = queries.get_all_signedoff_panels()

    td_data = td_data.apply(
        lambda row: checker.check_target(row, hgnc_data), axis=1
    )
    td_data = td_data.apply(
        lambda row: checker.check_test_method(row, td_config), axis=1
    )

    td_data = td_data.reindex(
        columns=[
            "Test ID", "Clinical Indication", "Target/Genes",
            "Identified targets", "Test Method", "Potential new test methods",
            "Potential removed test methods"
        ]
    )

    checker.compare_gp_td(td_data, genepanels_data, hgnc_data, signedoff_panels)

    td_data.to_excel("Checked_td.xlsx", sheet_name="data", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Checks the content of the rare disease sheet of the test "
            "directory"
        )
    )
    parser.add_argument(
        "test_directory", help="Test directory excel file to check"
    )
    parser.add_argument(
        "hgnc_dump", help="HGNC dump downloaded from the genenames.org website"
    )
    parser.add_argument(
        "genepanels",
        help="Genepanels file to compare against the provided test directory"
    )
    parser.add_argument(
        "-c", "--config", help="Test directory parser config file"
    )
    args = vars(parser.parse_args())
    main(args)
