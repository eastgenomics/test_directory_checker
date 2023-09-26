import numpy as np
import pandas as pd
from panelapp import queries

from test_directory_checker import checker


def test_compare_gp_td(td_data, genepanels_data, hgnc_dump):
    td_data = pd.read_csv(
        td_data, sep="\t", keep_default_na=False,
        converters={
            "Identified panels": lambda x: x.strip("[]").split(", ") if x != "" else [],
            "Identified genes": lambda x: x.strip("[]").split(", ") if x != "" else [],
        }
    )
    genepanels_data = pd.read_csv(
        genepanels_data, sep="\t", names=["ci", "panel", "gene"]
    )
    hgnc_dump = pd.read_csv(hgnc_dump, sep="\t")
    signedoff_panels = queries.get_all_signedoff_panels()

    expected_identical_cis = pd.DataFrame(
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

    expected_replaced_cis = pd.DataFrame(
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

    identical_cis, replaced_cis = checker.compare_gp_td(
        td_data, genepanels_data, hgnc_dump, signedoff_panels
    )

    for col in identical_cis.columns:
        np.testing.assert_array_equal(
            identical_cis[col].values, expected_identical_cis[col].values
        )

    for col in replaced_cis.columns:
        np.testing.assert_array_equal(
            replaced_cis[col].values, expected_replaced_cis[col].values
        )
