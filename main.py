import argparse

from test_directory_checker import checker, utils


def main(args):
    config = utils.load_config("configs/column_config.json")
    td_data = utils.parse_td(args["test_directory"], config)
    hgnc_data = utils.parse_hgnc_dump(args["hgnc_dump"])
    targets = checker.check_targets(
        td_data.to_dict()["Target/Genes"], hgnc_data
    )
    utils.write_json(targets)
    # checker.check_gene_symbols(td_data, hgnc_data)
    # checker.check_test_methods(td_data["Test Method"])



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
    args = vars(parser.parse_args())
    main(args)
