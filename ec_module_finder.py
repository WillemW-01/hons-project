import requests
import sys

from Set import Set

if __name__ == "__main__":
    ecs = Set()

    # fmt: off
    if len(sys.argv) == 1:
      ecs.add_items([
        "1.1.1.100", "1.1.1.205","1.1.1.25","1.1.1.44","1.1.1.49","1.1.1.88","1.1.1.94","1.2.4.4","1.3.1.10","1.3.1.98","2.1.1.185","2.1.1.45", 
        "2.1.3.2", "2.2.1.1", "2.3.1.266", "2.3.1.267", "2.3.1.39", "2.3.1.41", "2.3.1.51", "2.4.1.227", "2.4.2.10", "2.5.1.15", "2.5.1.19", 
        "2.5.1.31", "2.5.1.61", "2.5.1.7", "2.6.1.16", "2.7.1.33", "2.7.1.71","2.7.2.3","2.7.4.14","2.7.4.2","2.7.4.3","2.7.6.3",
        "2.7.6.5", "2.7.7.101", "2.7.7.23", "2.7.7.3", "2.7.7.39", "2.7.7.60", "2.7.8.13", "3.5.2.3", "3.5.3.1", "3.6.1.1", "3.6.5.3", "4.1.1.23",
        "4.1.1.33", "4.1.1.37", "4.1.1.71", "4.1.2.25", "4.1.3.36", "4.2.1.10", "4.2.1.24", "4.2.1.75", "4.2.3.4", "4.2.3.5", "4.98.1.1", "5.1.1.3",
        "5.1.3.1", "5.3.3.2", "5.4.2.11", "6.1.1.13", "6.2.1.5", "6.3.1.2", "6.3.2.13", "6.3.2.17", "6.3.2.4", "6.3.2.8",
        "6.3.2.9","6.3.4.15","6.3.4.2","6.3.5.2","6.3.5.5","6.4.1.2 "
      ])
    # fmt: on
    else:
        with open(sys.argv[1]) as ec_file:
            for line in ec_file:
                ecs.add(line.strip().lower())

    for ec in ecs:
        response = requests.get(f"https://rest.kegg.jp/get/{ec}")

        pathways = []
        should_add = False
        for line in response.text.split("\n"):
            if line.startswith("PATHWAY"):
                should_add = True

            if line.startswith("ORTHOLOGY"):
                should_add = False

            if should_add:
                pathways.append(line)

        print(f"Pathways for {ec}:")
        for pathway in pathways:
            print(f"{pathway}")
