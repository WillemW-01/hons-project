import requests
import re
import os
import time
import sys


def print_progress(i, length):
    completed = i / length
    not_completed = 80 - int(80 * completed)
    hash_part = "#" * int(80 * completed)
    empty_part = " " * (not_completed)

    print(f"[{hash_part}{empty_part}] {completed * 100:.2f}%", end="\r")
    time.sleep(0.01)


def write_set(uniq_set, output_file):
    for i, rxn in enumerate(uniq_set):
        if i < len(uniq_set):
            output_file.write(f'"{rxn}",')
        else:
            output_file.write(f'"{rxn}"],\n')


module_filepath = sys.argv[1]
output_path = sys.argv[2]

with open(module_filepath) as input_file:
    # First, extract the module ids from each line
    modules = []
    for line in input_file:
        matched = re.match(r"\s+(M\d{5})", line)
        if matched is not None:
            # the module id will be contained in match group 1
            modules.append(matched[1])

    # Then, extract all reaction ids for each module
    all_rxns = set()
    with open(output_path, "x") as output_file:
        output_file.write("{\n")
        for i, module in enumerate(modules):
            reactions = set()
            response = requests.get(f"https://rest.kegg.jp/get/{module}")
            if response.ok:
                text = response.text
                should_add = False
                temp_reactions = []
                for line in text.split("\n"):
                    if line.startswith("REACTION"):
                        should_add = True
                    if line.startswith("COMPOUND"):
                        should_add = False
                    if should_add:
                        matched = re.findall(r"R\d{5}", line)
                        if matched is not None:
                            temp_reactions += matched[0:]

                # add all reactions to a set to remove duplicates
                [reactions.add(reaction) for reaction in temp_reactions]
                [all_rxns.add(reaction) for reaction in temp_reactions]

                # able to keep track of progress in real time
                print_progress(i + 1, len(modules))
                output_file.write(f'"{module}":' + "[")
                if len(reactions) == 0:
                    output_file.write("],\n")
                for i, rxn in enumerate(reactions):
                    if i < len(reactions) - 1:
                        output_file.write(f'"{rxn}",')
                    else:
                        output_file.write(f'"{rxn}"],\n')
                output_file.flush()
                os.fsync(output_file.fileno())
            else:
                output_file.write(f"Error with {module}\n")

        # finally, print out the large, unique set of reactions
        output_file.write('"All": [')
        for i, rxn in enumerate(all_rxns):
            if i < len(all_rxns) - 1:
                output_file.write(f'"{rxn}",')
            else:
                output_file.write(f'"{rxn}"]\n')
        output_file.write("}\n")
