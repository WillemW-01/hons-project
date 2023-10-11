# Metabolic Network Writer

## Purpose

The code in this project serves to build a directed graph representing a metabolic network. This is done in order to visualize the differences between minimal metabolic networks of different organisms.

## Input

Running the file `main.py` requires the following paramers:

`$ python main.py <input_file> <job_name> [reaction_file]`

Where:

- `input_file` is a text file containing a list of EC numbers. These EC numbers can be in between other text as they are extracted using regex.
- `job_name` is appended to every file produced in the output.
- `reaction_file` is an optional file listing constraining the reactions that will be used for the output. These reactions must be KEGG REACTION identifiers (Rxxxxx).

## Output

The program will output the following:

- `output/job_name_nodes.gdf` the file that is used for visualization in the program Gephi.
- `output/job_name_edges.gdf` the file that represents the edges between nodes. Appended to the `nodes.gdf` file during execution.
- `logs/job_name.log`: a log file for more information on what happened during execution.

## Usage

1. First install the two depencies with `pip3 install -r requirements.txt`
1. Run `main.py` with the required input paramers. Files in the `input/` folder are all valid.
1. [Download](https://gephi.org/users/download/) and open the program Gephi.
1. Create a new project and import the nodes file from `output/job_name_nodes.gdf`
