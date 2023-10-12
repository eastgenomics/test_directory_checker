from pathlib import Path

from jinja2 import Environment, FileSystemLoader
import pandas as pd

ROOT_DIR = Path(__file__).absolute().parents[0]


def output_plot(
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

    output_folder.mkdir(exist_ok=True)
    environment = Environment(loader=FileSystemLoader(
        ROOT_DIR.joinpath("template")
    ))
    template = environment.get_template("table_template.html")
    content = template.render(
        filtered_tables=[table.to_html() for table in filtered_tables],
        title=output_name,
        full_table=full_table.to_html()
    )

    with open(
        f"{output_folder}/{output_name}", mode="w", encoding="utf-8"
    ) as f:
        f.write(content)
