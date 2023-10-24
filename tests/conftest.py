import pytest


def pytest_addoption(parser):
    parser.addoption("--td")
    parser.addoption("--genepanels")
    parser.addoption("--hgnc_dump")
    parser.addoption("--config")
    parser.addoption("--blacklist_config")


@pytest.fixture
def td_data(request):
    return request.config.getoption("--td")


@pytest.fixture
def genepanels_data(request):
    return request.config.getoption("--genepanels")


@pytest.fixture
def hgnc_dump(request):
    return request.config.getoption("--hgnc_dump")


@pytest.fixture
def config(request):
    return request.config.getoption("--config")


@pytest.fixture
def blacklist_config(request):
    return request.config.getoption("--blacklist_config")
