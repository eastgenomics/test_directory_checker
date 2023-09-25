import pytest


def pytest_addoption(parser):
    parser.addoption("--td")
    parser.addoption("--genepanels")
    parser.addoption("--hgnc_dump")


@pytest.fixture
def td_data(request):
    return request.config.getoption("--td")


@pytest.fixture
def genepanels_data(request):
    return request.config.getoption("--genepanels")


@pytest.fixture
def hgnc_dump(request):
    return request.config.getoption("--hgnc_dump")
