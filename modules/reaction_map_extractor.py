import json
import requests
import re
import time


def print_progress(i, length):
    completed = i / length
    not_completed = 80 - int(80 * completed)
    hash_part = "#" * int(80 * completed)
    empty_part = " " * (not_completed)

    print(f"[{hash_part}{empty_part}] {completed * 100:.2f}%", end="\r")
    time.sleep(0.01)


# fmt: off
reaction_maps = json.load(open("aureus_reaction_maps.json"))
total_len = sum([len(reaction_maps[item]) for item in reaction_maps])

with open("aureus_reactions_new.txt", "w") as output_file:
  i = 0
  for module in reaction_maps.keys():
      rxn_map_list = reaction_maps[module]
      output_file.write(f"> {module}\n")
      for rxn_map in rxn_map_list:
          response = requests.get(f"https://rest.kegg.jp/get/{rxn_map}")
          if response.ok:
              text = response.text
              should_add = False
              rxn_list = []
              for line in text.split("\n"):
                  if line.startswith("REACTION"):
                      should_add = True
                  if line.startswith("COMPOUND") or line.startswith("REL_PATHWAY"):
                      should_add = False
                  if should_add:
                      matched = re.findall(r"R\d{5}", line)
                      if matched is not None:
                          rxn_list += matched[0:]
              [output_file.write(rxn + "\n") for rxn in rxn_list]
          else:
              output_file.write(f"{rxn_map} not okay\n")
          i += 1
          print_progress(i, total_len)

# fmt: on
