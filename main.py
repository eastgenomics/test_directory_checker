import argparse
from pathlib import Path

from panelapp import queries

from test_directory_checker import checker, utils, output


def main(args):
    config = utils.load_config("configs/column_config.json")
    td_config = utils.load_config(args["config"])
    td_data = utils.parse_td(args["test_directory"], config)
    hgnc_data = utils.parse_hgnc_dump(args["hgnc_dump"])
    genepanels_data = utils.parse_genepanels(args["genepanels"])
    signedoff_panels = queries.get_all_signedoff_panels()

    target_data = td_data.apply(
        lambda row: checker.check_target(row, hgnc_data), axis=1
    )
    test_method_data = td_data.apply(
        lambda row: checker.check_test_method(row, td_config), axis=1
    )

    target_data = target_data.reindex(
        columns=[
            "Test ID", "Clinical Indication", "Target/Genes",
            "Identified panels", "Identified genes", "Test Method"
        ]
    )

    gene_locus_type = utils.get_locus_status_genes(
        target_data, signedoff_panels, hgnc_data
    )

    test_method_data = test_method_data.reindex(
        columns=[
            "Test ID", "Clinical Indication", "Test Method",
            "Potential new test methods"
        ]
    )

    (
        identical_tests, removed_tests, replaced_tests
    ) = checker.compare_gp_td(
        target_data, genepanels_data, hgnc_data, signedoff_panels, gene_locus_type
    )

    new_cis = checker.find_new_clinical_indications(
        target_data, genepanels_data
    )

    for df in [new_cis, target_data, test_method_data]:
        df.sort_values(["Test Method", "Test ID"], inplace=True)

    genes_to_check = utils.get_genes_from_td_target(
        target_data, signedoff_panels, gene_locus_type
    )

    presence_db_df = checker.check_if_genes_present_in_db(
        args["db_user"], args["db_password"], args["db_name"], "mysql",
        genes_to_check
    )

    # outputting logic
    output_folder = Path(args["output"])

    # filter tests have None in both the removed and added columns
    filtered_df = utils.filter_out_df(
        identical_tests, removed=None, added=None
    )
    output.output_table(
        identical_tests, "identical_tests.html", output_folder, filtered_df
    )

    output.output_table(removed_tests, "removed_tests.html", output_folder)

    # filter out tests have None in both the removed and added columns
    filtered_df = utils.filter_out_df(replaced_tests, removed=None, added=None)
    output.output_table(
        replaced_tests, "replaced_tests.html", output_folder, filtered_df
    )

    # filter out tests that have empty lists in the Identified panels and
    # Identified genes
    filtered_df = target_data.loc[
        (
            target_data["Identified panels"].str.len() == 0
        ) &
        (
            target_data["Identified genes"].str.len() == 0
        )
    ]
    output.output_table(
        target_data, "targets.html", output_folder, filtered_df
    )

    # filter out tests that have an empty string in the Potential new test
    # methods column
    filtered_df = utils.filter_out_df(
        test_method_data, **{"Potential new test methods": ""}
    )
    output.output_table(
        test_method_data, "test_methods.html", output_folder, filtered_df
    )

    # filter using tests that are present in the test methods in the passed
    # config file
    filtered_df = new_cis[
        new_cis["Test Method"].isin(td_config["ngs_test_methods"])
    ]
    output.output_table(new_cis, "new_cis.html", output_folder, filtered_df)

    # filter to get tests that have False in the presence_in_db or
    # has_clinical_transcript columns
    filtered_df = presence_db_df[
        (
            ~presence_db_df["presence_in_db"]
        ) |
        (
            ~presence_db_df["has_clinical_transcript"]
        )
    ]
    output.output_table(
        presence_db_df, "presence_in_db.html", output_folder, filtered_df
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Checks the content of the rare disease sheet of the test "
            "directory, more information in "
            "https://github.com/eastgenomics/test_directory_checker"
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

    parser.add_argument("db_user", help="Username of the database to check")
    parser.add_argument(
        "db_password", help="Username's password of the database to check"
    )
    parser.add_argument("db_name", help="Name of the database to check")

    parser.add_argument(
        "-c", "--config", help="Test directory parser config file"
    )
    parser.add_argument(
        "-o", "--output", help="Output folder", default="td_checker_output"
    )
    args = vars(parser.parse_args())
    main(args)
