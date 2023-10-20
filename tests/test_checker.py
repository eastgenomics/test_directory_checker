import json
from pathlib import Path

import numpy as np
import pandas as pd
from panelapp import queries
import pytest

from test_directory_checker import checker, output, utils


@pytest.fixture
def setup_hgnc_dump(hgnc_dump):
    yield pd.read_csv(hgnc_dump, sep="\t")


@pytest.fixture
def setup_config(config):
    config = open(config)
    data = json.load(config)
    config.close()
    yield data


@pytest.fixture
def setup_td_data(td_data):
    yield pd.read_csv(
        td_data, sep="\t", keep_default_na=False,
        converters={
            "Identified panels": lambda x: x.strip("[]").split(", ") if x != "" else [],
            "Identified genes": lambda x: x.strip("[]").split(", ") if x != "" else [],
        }
    )


@pytest.fixture
def setup_genepanels_data(genepanels_data):
    yield pd.read_csv(
        genepanels_data, sep="\t", names=["ci", "panel", "gene"]
    )


@pytest.fixture
def setup_signedoff_panels():
    yield queries.get_all_signedoff_panels()


def test_check_target_panel(setup_hgnc_dump):
    """ Test for checking a panel target i.e. Foo (123) should find 123 as a
    panel target

    Args:
        setup_hgnc_dump (function): Fixture that parsed the hgnc dump
    """

    row = pd.Series(
        ["Thoracic aortic aneurysm or dissection (700)"],
        index=["Target/Genes"]
    )
    processed_row = checker.check_target(row, setup_hgnc_dump)
    expected_row = pd.Series(
        [
            "Thoracic aortic aneurysm or dissection (700)",
            ["700"],
            []
        ],
        index=["Target/Genes", "Identified panels", "Identified genes"]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_check_target_gene(setup_hgnc_dump):
    """ Test for checking a panel target i.e. BRCA1, BRCA2 should return a list
    containing HGNC:1100 and HGNC:1101

    Args:
        setup_hgnc_dump (function): Fixture that parses the hgnc dump
    """

    row = pd.Series(
        ["BMPR2"],
        index=["Target/Genes"]
    )
    processed_row = checker.check_target(row, setup_hgnc_dump)
    expected_row = pd.Series(
        [
            "BMPR2",
            [],
            ["HGNC:1078"]
        ],
        index=["Target/Genes", "Identified panels", "Identified genes"]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_check_test_method_exists(setup_config):
    """ Test that finds an existing test method

    Args:
        setup_config (function): Fixture that loads a JSON config file
    """

    row = pd.Series(["Small panel"], index=["Test Method"])
    processed_row = checker.check_test_method(row, setup_config)
    expected_row = pd.Series(
        [
            "Small panel", ""
        ], index=[
            "Test Method", "Potential new test methods"
        ]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_check_test_method_not_exists(setup_config):
    """ Test that finds a new test method

    Args:
        setup_config (function): Fixture that loads a JSON config file
    """

    row = pd.Series(["New test method"], index=["Test Method"])
    processed_row = checker.check_test_method(row, setup_config)
    expected_row = pd.Series(
        [
            "New test method", "New test method"
        ], index=[
            "Test Method", "Potential new test methods",
        ]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_compare_gp_td(
    setup_td_data, setup_genepanels_data, setup_hgnc_dump,
    setup_signedoff_panels
):
    """ Test to check the output of the compare_gp_td function. 3 manually
    written dataframes check the 3 outputs of the function

    Args:
        setup_td_data (function): Fixture that parses the test directory data
        setup_genepanels_data (function): Fixture that parses the genepanels data
        setup_hgnc_dump (function): Fixture that parses the hgnc dump
        setup_signedoff_panels (function): Fixture that creates the signedoff
        panel dictionary
    """

    gene_locus_type = utils.get_locus_status_genes(
        setup_td_data, setup_signedoff_panels, setup_hgnc_dump
    )

    expected_identical_tests = pd.DataFrame(
        [
            [
                "R130.1_Short QT syndrome_P", "Short QT syndrome_3.1",
                "HGNC:1390, HGNC:6251, HGNC:6263, HGNC:6294", "R130.1",
                "Short QT syndrome (224)", "3.1",
                "HGNC:1390, HGNC:6251, HGNC:6263, HGNC:6294",
                None, None
            ],
            [
                "R341.1_Hereditary angioedema types I and II_G",
                "HGNC:1228_SG_panel_1.0.0", "HGNC:1228", "R341.1",
                "SERPING1", "",
                "HGNC:1228",
                None, None
            ],
            [
                "R122.1_Factor XIII deficiency_P",
                "HGNC:3531_SG_panel_1.0.0, HGNC:3534_SG_panel_1.0.0",
                "HGNC:3531, HGNC:3534", "R122.1",
                "F13A1, F13B, F2", "",
                "HGNC:3531, HGNC:3534, HGNC:6666",
                None, "HGNC:6666"
            ],
            [
                "R143.1_Neonatal diabetes_P",
                "HGNC:59_SG_panel_1.0.0, HGNC:6257_SG_panel_1.0.0",
                "HGNC:59, HGNC:6257", "R143.1",
                "ABCC8", "",
                "HGNC:59",
                "HGNC:6257", None
            ],
            [
                "R347.1_Inherited predisposition to acute myeloid leukaemia (AML)_P",
                "Inherited predisposition to acute myeloid leukaemia (AML)_3.0",
                "HGNC:10471, HGNC:11730, HGNC:11998, HGNC:1349, HGNC:1833, HGNC:18674, HGNC:29186, HGNC:3495, HGNC:4171, HGNC:6919",
                "R347.1",
                "Inherited predisposition to acute myeloid leukaemia (AML) (525)",
                "3.0",
                "HGNC:10471, HGNC:11730, HGNC:11998, HGNC:1349, HGNC:1833, HGNC:18674, HGNC:29186, HGNC:3495, HGNC:4171, HGNC:6919",
                None, None
            ],
        ],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    expected_removed_tests = pd.DataFrame(
        [
            [
                "R1000.1_Removed test", "Removed test panel",
                "HGNC:1390, HGNC:6251, HGNC:6263, HGNC:6294"
            ]
        ],
        columns=[
            "gemini_name", "panel", "genes"
        ]
    )

    expected_replaced_tests = pd.DataFrame(
        [
            [
                "R100.1_Test clinical indication", "Test panel_1.0",
                "HGNC:10471, HGNC:11730, HGNC:11998, HGNC:1349, HGNC:1833, HGNC:18674, HGNC:29090, HGNC:29186, HGNC:3495, HGNC:4171, HGNC:6919",
                "R100.2", "Test panel",
                "3.0",
                "HGNC:10471, HGNC:11730, HGNC:11998, HGNC:1349, HGNC:1833, HGNC:18674, HGNC:29186, HGNC:3495, HGNC:4171, HGNC:6919",
                "HGNC:29090", None
            ],
            [
                "R134.1_Familial hypercholesterolaemia_P",
                "Familial hypercholesterolaemia (GMS)_2.0",
                "HGNC:20001, HGNC:603, HGNC:613, HGNC:6547",
                "R134.2", "Familial hypercholesterolaemia", "2.0",
                "HGNC:18640, HGNC:20001, HGNC:603, HGNC:613, HGNC:6547",
                None, "HGNC:18640"
            ],
            [
                "R134.1_Familial hypercholesterolaemia_P",
                "Familial hypercholesterolaemia (GMS)_2.0",
                "HGNC:20001, HGNC:603, HGNC:613, HGNC:6547",
                "R134.3", "Familial hypercholesterolaemia", "",
                "HGNC:1228",
                "HGNC:20001, HGNC:603, HGNC:613, HGNC:6547", "HGNC:1228"
            ]
        ],
        columns=[
            "gemini_name", "panel", "genes", "td_ci", "td_target",
            "td_version", "td_genes", "removed", "added"
        ]
    )

    (
        identical_tests, removed_tests, replaced_tests
    ) = checker.compare_gp_td(
        setup_td_data, setup_genepanels_data, setup_signedoff_panels,
        gene_locus_type
    )

    for col in identical_tests.columns:
        np.testing.assert_array_equal(
            identical_tests[col].to_numpy(),
            expected_identical_tests[col].to_numpy()
        )

    for col in removed_tests.columns:
        np.testing.assert_array_equal(
            removed_tests[col].to_numpy(),
            expected_removed_tests[col].to_numpy()
        )

    for col in replaced_tests.columns:
        np.testing.assert_array_equal(
            replaced_tests[col].to_numpy(),
            expected_replaced_tests[col].to_numpy()
        )

    # not part of the actual test but to provide examples in the readme
    # screenshots
    filtered_df = utils.filter_out_df(
        identical_tests, removed=None, added=None
    )
    output.output_table(
        identical_tests, "identical_tests.html", Path("tests/test_outputs"),
        filtered_df
    )

    output.output_table(
        removed_tests, "removed_tests.html", Path("tests/test_outputs")
    )

    filtered_df = utils.filter_out_df(replaced_tests, removed=None, added=None)
    output.output_table(
        replaced_tests, "replaced_tests.html", Path("tests/test_outputs"),
        filtered_df
    )


def test_find_new_clinical_indications(setup_td_data, setup_genepanels_data):
    """ Test to find new clinical indications. 4 bespoke tests were added to
    the test directory data that the code is supposed to pick up

    Args:
        setup_td_data (function): Fixture that parses the test directory data
        setup_genepanels_data (function): Fixture that parses the genepanels data
    """

    expected_new_cis = pd.DataFrame(
        [
            [
                "R100.2", "Test clinical indication", "Test panel", ["525"],
                [], "Small panel", "set()"
            ],
            [
                "R134.2", "Familial hypercholesterolaemia",
                "Familial hypercholesterolaemia", ["772"], [],
                "Small panel", "set()"
            ],
            [
                "R134.3", "Familial hypercholesterolaemia",
                "Familial hypercholesterolaemia", [], ["HGNC:1228"],
                "Single gene sequencing <=10 amplicons", "set()"
            ],
            [
                "R500.1", "Test new ci", "Test new panel", ["525"],	[],
                "Small panel", "set()"
            ],
            [
                "R666.1", "Test no clinical transcript",
                "Test no clinical transcript", [], ["HGNC:9999"],
                "Single gene sequencing <=10 amplicons", "set()"
            ]
        ],
        columns=[
            "Test ID", "Clinical Indication", "Target/Genes",
            "Identified panels", "Identified genes", "Test Method",
            "Potential new test methods"
        ]
    )

    new_cis = checker.find_new_clinical_indications(
        setup_td_data, setup_genepanels_data
    )

    for col in new_cis.columns:
        np.testing.assert_array_equal(
            new_cis[col].values, expected_new_cis[col].values
        )


def test_check_if_genes_in_db(
    setup_td_data, setup_hgnc_dump, setup_signedoff_panels
):
    gene_locus_type = utils.get_locus_status_genes(
        setup_td_data, setup_signedoff_panels, setup_hgnc_dump
    )

    genes_to_check = utils.get_genes_from_td_target(
        setup_td_data, setup_signedoff_panels, gene_locus_type
    )

    presence_in_db_df = checker.check_if_genes_present_in_db(
        "", "", "tests/test_files/test_db.db", "sqlite", genes_to_check
    )

    expected_presence_genes = pd.DataFrame(
        [
            ["HGNC:10471", True, True],
            ["HGNC:11730", True, True],
            ["HGNC:11998", True, True],
            ["HGNC:1228", True, True],
            ["HGNC:1349", True, True],
            ["HGNC:1390", True, True],
            ["HGNC:1833", True, True],
            ["HGNC:18640", True, True],
            ["HGNC:18674", True, True],
            ["HGNC:20001", True, True],
            ["HGNC:29186", True, True],
            ["HGNC:3495", True, True],
            ["HGNC:3531", True, True],
            ["HGNC:3534", True, True],
            ["HGNC:4171", True, True],
            ["HGNC:59", True, True],
            ["HGNC:603", True, True],
            ["HGNC:613", True, True],
            ["HGNC:6251", True, True],
            ["HGNC:6263", True, True],
            ["HGNC:6294", True, True],
            ["HGNC:6547", True, True],
            ["HGNC:6666", False, False],
            ["HGNC:6919", True, True],
            ["HGNC:9999", True, False]
        ],
        columns=["gene", "presence_in_db", "has_clinical_transcript"]
    )

    for col in presence_in_db_df.columns:
        np.testing.assert_array_equal(
            presence_in_db_df[col].to_numpy(),
            expected_presence_genes[col].to_numpy()
        )
