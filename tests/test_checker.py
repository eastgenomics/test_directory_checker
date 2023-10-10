import json

import numpy as np
import pandas as pd
from panelapp import queries
import pytest

from test_directory_checker import checker


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
    row = pd.Series(["Small panel"], index=["Test Method"])
    processed_row = checker.check_test_method(row, setup_config)
    expected_row = pd.Series(
        [
            "Small panel", set()
        ], index=[
            "Test Method", "Potential new test methods"
        ]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_check_test_method_not_exists(setup_config):
    row = pd.Series(["New test method"], index=["Test Method"])
    processed_row = checker.check_test_method(row, setup_config)
    expected_row = pd.Series(
        [
            "New test method", {"New test method"}
        ], index=[
            "Test Method", "Potential new test methods",
        ]
    )

    np.testing.assert_array_equal(processed_row, expected_row)


def test_compare_gp_td(
    setup_td_data, setup_genepanels_data, setup_hgnc_dump,
    setup_signedoff_panels
):
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
                "R344.1_Primary hyperaldosteronism - KCNJ5_G",
                "HGNC:6266_SG_panel_1.0.0", "HGNC:6266", "R344.1",
                "KCNJ5", "",
                "HGNC:6266",
                None, None
            ],
            [
                "R345.2_Facioscapulohumeral muscular dystrophy - extended testing_G",
                "HGNC:29090_SG_panel_1.0.0", "HGNC:29090", "R345.2",
                "SMCHD1", "",
                "HGNC:29090",
                None, None
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
            "gemini_name", "panel", "genes", "td_ci", "td_target", "td_version",
            "td_genes", "removed", "added"
        ]
    )

    identical_tests, removed_tests, replaced_tests = checker.compare_gp_td(
        setup_td_data, setup_genepanels_data, setup_hgnc_dump,
        setup_signedoff_panels
    )

    for col in identical_tests.columns:
        np.testing.assert_array_equal(
            identical_tests[col].values, expected_identical_tests[col].values
        )

    for col in removed_tests.columns:
        np.testing.assert_array_equal(
            removed_tests[col].values, expected_removed_tests[col].values
        )

    for col in replaced_tests.columns:
        np.testing.assert_array_equal(
            replaced_tests[col].values, expected_replaced_tests[col].values
        )


def test_find_new_clinical_indications(setup_td_data, setup_genepanels_data):
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
                "R500.1", "Test new ci", "Test new panel", ["100"],	[],
                "Small panel", "set()"
            ],
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
