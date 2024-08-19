# Antibody Discovery and Ortholog Identification Script

## Introduction

This script is designed to facilitate the discovery of antibodies and the identification of orthologs across different species. The script aligns protein sequences, maps gene IDs, and searches for ortholog genes in target species. It processes phosphorylation site datasets to provide detailed outputs, including aligned protein sequences, ortholog information, and analysis of phosphorylation sites.

This tool is particularly useful for researchers working in the fields of bioinformatics, genomics, and proteomics, where cross-species comparison of gene and protein data is essential.

## Requirements

The script requires the following dependencies:

- Python 3.7+
- pandas
- Biopython
- Openpyxl
- Additional custom modules:
  - `Logger`
  - `TimeTracker`

These dependencies can be installed via `pip` if they are not already available in your Python environment.

## Installation

To set up the environment and install the necessary dependencies, follow these steps:

1. **Clone the repository:**

    ```bash
    git clone git@github.com:chenxi-gao/antibody_discovery.git
    cd antibody-discovery-script
    ```

2. **Create and activate a virtual environment (optional but recommended):**

    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. **Install the required Python packages:**

    ```bash
    pip install -r requirements.txt
    ```

4. **Ensure that the `Logger` and `TimeTracker` modules are available in your Python path. If they are not included in the repository, add them or install them, according to the provided documentation.**

## Usage

### Running the Script

To run the script, use the following command:

```bash
python antibody_discovery.py [species] [input_file_path]
