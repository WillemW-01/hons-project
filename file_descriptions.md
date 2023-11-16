# Descriptions of files in this repository

## Folder descriptions

| Folder          | Description                                                                                                                    |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `/cached_ecs`   | Hosts text for EC numbers so that requests need not be sent for commonly found enzymes (in the datasets).                      |
| `/calculations` | Holds the excel files and other files used for module comparison and minimal set identification.                               |
| `/gephi`        | Projects in the Gephi program are saved in this folder. One can open these in Gephi to view my own projects used for analysis. |
| `/images`       | The images used in the report and seminar. Generated using Gephi.                                                              |
| `/input`        | The sample input files used for generating the networks used in analysis.                                                      |
| `/logs`         | Holds a record of what the program did during execution. Useful for debugging.                                                 |
| `/modules`      | Groups files that are used for module assignment and some of the reaction data cache.                                          |
| `/output`       | The `.gdf` files are outputted to this directory. Use these in Gephi.                                                          |
| `/wcm`          | A last-minute test folder for whole-cell model data using the BiGG database. Only for those really interested.                 |

## `/.`

_(Top level folder)_

| File                   | Description                                                                                                                              |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `.gitignore`           | Used for ignoring certain files when using version control.                                                                              |
| `file_descriptions.md` | The current file. Explains the meaning behind file naming and uses of certain folders and files.                                         |
| `LICENSE`              | The file describing the limitations of distributing and reusing the code written in this repository. MIT license - very liberal license. |
| `main.py`              | The main utility used to execute the function of the program. Usage described in README                                                  |
| `README.md`            | File describing general information about the program, its installation and use.                                                         |
| `requirements.txt`     | This file is required to install the necessary python packages for normal functioning. See README.                                       |
| `Set.py`               | A wrapper class for a custom implementation of the built-in `set` object.                                                                |

## `/calculations`

| File                            | Description                                                                                                                                                 |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Comparison.xlsx`               | This spreadsheet hosts the majority of the comparisons made for the minimal set identification. Sheets exist for both enzyme and reaction comparisons.      |
| `Modularity Statistics.xlsx`    | This spreadsheet contains the data proving that after removing energy molecules from the graphs, that the graphs indeed became less connected.              |
| `/comparison/*_ecs.txt`         | These files contains all the ec numbers of the enzymes in the final graphs, they were used to do set operations in order to find intersections and unions.  |
| `/comparison/*_rxns.csv`        | These files contain reaction identifiers that were used for the same purpose as mentioned above, except for comparing reactions.                            |
| `/comparison/set_operations.py` | This file was used to perform the set operations (intersections and unions) to compare the different datasets, using the `*_ecs.txt` and `*_rnxs.csv` files |

## `/gephi`

| File          | Description                                                           |
| ------------- | --------------------------------------------------------------------- |
| `final.gephi` | This file was used to generate the figures in the report and seminar. |
| `full.gephi`  | This file holds the experimental WCM approach's graphs                |

## `/images`

| File          | Description                                                                              |
| ------------- | ---------------------------------------------------------------------------------------- |
| `/*_after.*`  | These images are the versions of the graphs after pruning of energy molecules were done  |
| `/*_before.*` | These images are the versions of the graphs before pruning of energy molecules were done |

## `/input`

| File                             | Description                                                                                                                        |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `full.tsv`                       | Contains the full genome and essentiality annotation of the Agrobacterium Tumafaciens organism. Used in the beginning for testing. |
| `jcvi_ecs.txt`                   | Contains a curated list of EC numbers for the JCVI dataset.                                                                        |
| `jcvi_reactions.txt`             | Contains a curated list of KEGG reaction identifiers for the JCVI dataset.                                                         |
| `M_genitalium_ecs_essential.txt` | Contains the ec numbers marked as essential for M genitalium                                                                       |
| `M_genitalium_ecs_full.txt`      | Contains all the ec numbers for M genitalium                                                                                       |
| `M_pneumoniae_ecs_essential.txt` | Contains the ec numbers marked as essential for M pneumoniae                                                                       |
| `M_pneumoniae_ecs_full.txt`      | Contains all the ec numbers for M pneumoniae                                                                                       |
| `Minesweeper_ecs_essential.txt`  | Contains the ec numbers marked as essential for Minesweeper (M genitalium based)                                                   |
| `Minesweeper_ecs_full.txt`       | Contains all the ec numbers for Minesweeper (M genitalium based)                                                                   |
| `S_aureus_ecs_essential.txt`     | Contains the ec numbers marked as essential for S Aureus                                                                           |
| `S_aureus_ecs_full.txt`          | Contains all the ec numbers for S Aureus                                                                                           |
| `theoretical.txt`                | A test dataset that contains complete versions of the glycolysis, PPP, phospholipid and nucleotide metabolism modules.             |

## `/modules`

| File                            | Description                                                                                                                                             |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `all_modules.json`              | Unused file that was originally meant for assigning modules to enzymes in a previous version.                                                           |
| `all_reactions.json`            | This file holds the reaction maps that's used to assign modules to enzymes via the reactions they catalyse.                                             |
| `aureus_modules.json`           | This file was originally used to attempt and create a restricted list of reactions when generating S aureus output. It was not successful, however.     |
| `aureus_reaction_maps.json`     | Same as above.                                                                                                                                          |
| `aureus_reactions.json`         | Same as above.                                                                                                                                          |
| `cached_reaction_ids.json`      | Saves a list of EC numbers with the set of reactions that they catalyse. Used to speed up execution of the program.                                     |
| `cached_reactions.json`         | Saves a union of all the reactions in all the datasets. Used to speed up execution of the program.                                                      |
| `cached_reaction_ids.json`      | Saves a list of EC numbers with the set of reactions that they catalyse. Used to speed up execution of the program.                                     |
| `genitalium_modules.json`       | This file was originally used to attempt and create a restricted list of reactions when generating M genitalium output. It was not successful, however. |
| `genitalium_reaction_maps.json` | Same as above                                                                                                                                           |
| `map_conversion.json`           | Originally used to assign modules to enzymes. Not used anymore.                                                                                         |
| `module_colors.json`            | File that saves the color assignment of the modules.                                                                                                    |
| `reaction_maps.json`            | Used to fetch the module an enzyme belongs to after its reaction map was determined.                                                                    |
