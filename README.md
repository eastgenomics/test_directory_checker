# test_directory_checker

Suite of tests to check the content of the test directory.

## How to setup

### Setup the environment

```bash
python -m venv test_directory_checker
source test_directory_checker/bin/activate
pip install -r requirements.txt 
```

### Download the genepanels file

Find the genepanels that you want to use in DNAnexus in the `001_Reference` project in the `dynamic_files/gene_panels` folder.

### Download the test directory

Either get it from the bioinformatics team or download it from https://www.england.nhs.uk/publication/national-genomic-test-directories/
The version from the website might not be the latest one.

### Download the HGNC dump

From https://www.genenames.org/download/custom/, check the boxes as shown in the screenshot below.

![Alt text](image.png)

### Find the test directory config

The config file needed for the test directory checker is located here: https://github.com/eastgenomics/test_directory_parser/tree/main/configs

## How to run

```bash
source test_directory_checker/bin/activate
python3 main.py ${excel_td} ${hgnc_dump} ${genepanels_file} -c ${td_config}
```

## Tests

### Check the targets

Using the test directory, the code gets the targets and extracts the list of genes.
It outputs an HTML representation of the dataframe containing the identified genes/panels.

Output: `targets.html`

![Alt text](image-1.png)

### Check the test methods

Using the test directory, the code compares the test methods in the test directory and the ones in our config file that mimics the ones in the test directory parser.
It outputs an HTML representation of the dataframe containing the comparison of the test methods.

Output: `test_methods.html`

![Alt text](image-2.png)

### Compare genepanels and test directory

Using the test directory and a genepanels file, the code compares the content of each clinical indication in the genepanels file and finds the corresponding test code or clinical indication code to find if the genepanels clinical indication contains the genepanels ci.

Output: `identical_tests.html` + `removed_tests.html` + `replaced_tests.html`

![Alt text](image-3.png)

### Find new clinical indications

The code will go through every test in the test directory and it will attempt to find tests missing from the genepanels file.

Output: `new_cis.html`

![Alt text](image-4.png)

## Unittesting

The unit tests check the functions of the `checker.py` file. These are the main functions used to check the content of the test directory.

```bash
py.test -s tests/test_checker.py --td tests/test_files/test_td.tsv --genepanels tests/test_files/test_genepanels.tsv --hgnc_dump ${hgnc_dump} --config tests/test_files/test_config.json
```
