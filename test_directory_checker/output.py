from pathlib import Path

from jinja2 import Environment, FileSystemLoader
import pandas as pd

from test_directory_checker.utils import check_if_output_folder_exists, get_date

ROOT_DIR = Path(__file__).absolute().parents[0]


def mkdir_output_folder(output_folder: Path):
    """ Create the output folder

    Args:
        output_folder (Path): Path to the output folder

    Returns:
        Path: Created output folder path
    """

    date = get_date()
    counter = 1

    # add to the counter until the folder path doesn't exist
    while check_if_output_folder_exists(output_folder, date, counter):
        counter += 1 

    folder_path = output_folder / f"{date}-{counter}"
    folder_path.mkdir()
    return folder_path


def log_command_line(output_folder: Path, command_line: str):
    """ Write a txt file containing the command line used to run the script

    Args:
        output_folder (Path): Path to the output folder
        command_line (str): Command line at runtime
    """

    with open(output_folder / "command_line.txt", "w") as f:
        f.write(command_line)


def output_test_methods(table: pd.Series, output_name: str, output_folder: Path):
    """ Output the test methods as a 2 column dataframe which will contain the
    test methods not present in the config and the associated test IDs

    Args:
        table (pd.Series): Series containing the test methods and the test IDs associated
        output_name (str): Output name of the file
        output_folder (Path): Output folder
    """

    environment = Environment(loader=FileSystemLoader(
        ROOT_DIR.joinpath("template")
    ))
    template = environment.get_template("table_template.html")
    content = template.render(
        filtered_tables=[table.to_frame().to_html()],
        title=output_name,
    )

    with open(output_folder / output_name, mode="w", encoding="utf-8") as f:
        f.write(content)


def output_table(
    full_table: pd.DataFrame, output_name: str, output_folder: Path,
    *filtered_tables
):
    """ Output HTML of the data

    Args:
        full_table (pd.DataFrame): Dataframe containing all the data for that
        output
        output_name (str): Name of the HTML page
        output_folder (Path): Folder of the output
        *filtered_tables: List of filtered down tables originating from
        full_table
    """

    environment = Environment(loader=FileSystemLoader(
        ROOT_DIR.joinpath("template")
    ))
    template = environment.get_template("table_template.html")
    content = template.render(
        filtered_tables=[table.to_html() for table in filtered_tables],
        title=output_name,
        full_table=full_table.to_html()
    )

    with open(output_folder / output_name, mode="w", encoding="utf-8") as f:
        f.write(content)
