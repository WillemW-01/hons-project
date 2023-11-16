# Metabolic Network Writer

## Purpose

The code in this project serves to build a directed graph representing a metabolic network. This is done in order to visualize the differences between minimal metabolic networks of different organisms.

## Background

In order to study minimal metabolic networks, one can compare the metabolic networks generated from minimal genomes. This program was written specifically for this purpose. The included input contains minimal genomes of a few organisms. The files contained in the `output` folder were the previously generated graph files of these minimal genomes. The file `final.gephi` in the `gephi` folder can also be opened in Gephi to view the final graph project. A list of EC numbers for any organism can be used as input.

Visual exploration of metabolic networks allows comparison and validation of these networks.

The graph structure here is a bipartite graph, where enzymes and metabolites are represented by nodes, and their association as edges (metabolites catalysed in the same reaction by the same enzyme will be connected to that enzyme).

"Branches" of nodes stemming from the main part of a graph indicates successive metabolic steps being taken.

The enzymes are colored according to the main function of the reactions they catalyse. This can be viewed in the "Appearance" tab in Gephi.

Gephi can also be used to calculate the modularity, diameter, radius, shortest path, and other graph statistics using the "Statistics" tab.

## Environment

- The script requires a ` python` environment above version `3.8`.
- The Gephi program uses a `Java` environment, but this comes packaged with the download.

## Input

Running the file `main.py` requires the following paramers:

`$ python main.py <input_file> <job_name> [reaction_file]`

Where:

- `input_file` is a text file containing a list of EC numbers. These EC numbers can be in between other text as they are extracted using regex. Note that partial EC numbers will be ignored.
- `job_name` is appended to every file produced in the output.
- `reaction_file` is an optional file listing constraining the reactions that will be used for the output. These reactions must be KEGG REACTION identifiers (Rxxxxx). It is also expected that these reactions appear in the reactions catalysed by the enzymes listed in the input file.

## Output

The program will output the following:

- `output/<job_name>_nodes.gdf` the file that is used for visualization in the program Gephi.
- `logs/<job_name>.log`: a log file for more information on what happened during execution.

## Usage

### Python Script

1. First install the depencies with `pip3 install -r requirements.txt`
1. Run `main.py` with the required input paramers. Files in the `input/` folder are all valid.
1. The program detects the number of unique EC numbers in the input. It will then want to break this up into "batches" such that it won't be rate-limited by the KEGG database. The recommended batch size is between 50-75. If there are more than 2 batches, the program will wait 2 minutes between batches in order to prevent being rate-limited.
1. The progress of the program will be visible with different stages. "Retrieval" refers to gathering all the enzyme entries from either the cached files or the database, "Extracting" refers to extracing substrates and products from all the reactions, "Writing" refers to where the final file is being written. Note that the log file can be viewed at any time, as the output is flushed as the program executes.
1. At the end, the time of execution is displayed, and the `output/<job_name>_nodes.gdf` and `logs/<job_name>.log` files are generated.
1. You can move on to exploration in Gephi.

### Gephi Exploration

1. First, [download](https://gephi.org/) and install Gephi.
1. After installation, open Gephi.
1. Then, in the file tab, choose "New Project" (alternatively, if a popup occurs, choose "New Project")
1. You will see a Workspace tab open up. Then click on "Data Laboratory" near the top.
1. Then choose the option to import a spreadsheet.
1. Choose the file you just generated in the previous section (or any of the provided outputs). Note it must be a `.gdf` file.
1. With the import wizard popup, make sure that the "Edge merge strategy" is "Average", and that you choose the option "New workspace". Click "Ok".
1. You will see the nodes and edges are loaded into the program. Now switch to the "Overview" tab.
1. You will now see the network in a random configuration. Click on the "Layout" tab to change this.
1. Choose the "Random layout" option and enter a size of "1". Click the green run button.
1. Now, switch to the "ForceAtlas 2" spatialization algorithm. Change the following parameters:

   - "Scaling": 10.0,
   - "Gravity": 10.0

1. Then click on run and wait for the network to stabilise. You can also save these parameters in a preset with a chosen name (bottom of layout screen).
1. After the network stabilised, click on the up icon on the bottom-right of the screen. This will open label configuration.
1. Switch to the "Labels" tab and check the following options: "Hide non-selected" and then "Node".
1. Now, it will be possible to view the labels of the nodes.
1. Use the following controls with the mouse (hard with trackpad):

   - Zoom: mousewheel
   - Panning: right mouse
   - Move nodes: left mouse

In order to get the same output as the images in the report, you must remove all energy molecules, coenzymes and cofactors. Examples like these include ATP, ADP, CoA, diphosphate, NADP, NADH, etc. Use your discretion. Also remove "orthophosphate". When you run the spatialization algorithm again, the network topology will be much clearer.
