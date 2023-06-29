import argparse

from test_directory_checker import checker, utils


def main(args):
    config = utils.load_config("configs/column_config.json")
    td_config = utils.load_config(args["config"])
    td_data = utils.parse_td(args["test_directory"], config)
    hgnc_data = utils.parse_hgnc_dump(args["hgnc_dump"])
    targets = checker.check_targets(
        td_data.to_dict()["Target/Genes"], hgnc_data
    )
    new_test_methods, removed_test_methods = checker.check_test_methods(
        td_data.to_dict()["Test Method"], td_config
    )
    utils.write_json(targets)
    utils.write_test_methods(new_test_methods, removed_test_methods)


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
        "-c", "--config", help="Test directory parser config file"
    )
    args = vars(parser.parse_args())
    main(args)
