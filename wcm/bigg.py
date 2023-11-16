import requests
import json
import sys
import time

input_genes = sys.argv[1]
input_model = sys.argv[2]
job_name = sys.argv[3]

INPUT_TEXT = open(f"{input_genes}").read().split("\n")
INPUT_TEXT = [row.split(",") for row in INPUT_TEXT]
MODEL_DATA = json.load(open(input_model))
response_file = open(f"{job_name}_response.json", "w")
cache_file = json.load(open(f"{job_name}_read_response.json"))
NO_RESPONSE_GENES = open(f"{job_name}_skip_genes.txt").read().splitlines()
FILTER_LIST = open(f"filter_list.txt").read().splitlines()
log_file = open(f"{job_name}.log", "w")

MET_START = 2000


def print_progress(i, length, delay=None):
    completed = i / length
    not_completed = 80 - int(80 * completed)
    hash_part = "#" * int(80 * completed)
    empty_part = " " * (not_completed)

    print(f"[{hash_part}{empty_part}] {completed * 100:.2f}%", end="\r")

    if delay is not None:
        time.sleep(delay)


def is_reversible(reaction_dict):
    lower = float(reaction_dict["lower_bound"])
    upper = float(reaction_dict["upper_bound"])

    if lower < 0.0 and upper > 0.0:
        return True
    elif lower == 0.0 and upper > 0.0:
        return False
    elif lower == 0.0 and upper == 0.0:
        return None


def extract_metabolites(metabolite_dict):
    substrates = []
    products = []

    for met in metabolite_dict:
        if metabolite_dict[met] > 0:
            products.append(met)
        else:
            substrates.append(met)

    return substrates, products


def append_reaction_data(total_rxns, response_rxns):
    for rxn in response_rxns:
        bigg_id = rxn["bigg_id"]
        model_rxn = MODEL_DATA["reactions"][bigg_id]
        reversible = is_reversible(model_rxn)
        substrates, products = extract_metabolites(model_rxn["metabolites"])
        total_rxns.append(
            {
                "id": bigg_id,
                "is_reversible": reversible,
                "substrates": substrates,
                "products": products,
            }
        )


def add_group(bag: set, items):
    for item in items:
        bag.add(item)


def get_index(item: str, bag: set):
    try:
        pos = list(bag).index(item)
        return pos, False
    except:
        bag.add(item)
        return len(bag) - 1, True


def filter_metabolites(substrates, products, filter_list):
    substrates = [sub for sub in substrates if sub not in filter_list]
    products = [prod for prod in products if prod not in filter_list]
    return substrates, products


def write_met(pos_1, pos_2, file_obj, is_reverse=False):
    file_obj.write(f"{pos_1:04d},{pos_2:04d},1.0,true\n")
    if is_reverse:
        file_obj.write(f"{pos_2:04d},{pos_1:04d},1.0,true\n")


def log(message):
    log_file.write(message + "\n")


def write_reaction(
    index: int, reaction: dict, metabolites: set, node_file, edge_file, met_file
):
    rxn_id = reaction["id"]
    substrates, products = reaction["substrates"], reaction["products"]
    substrates, products = filter_metabolites(substrates, products, FILTER_LIST)

    reversible = reaction["is_reversible"]
    if reversible != "None":
        node_file.write(f"{index:04d},{rxn_id},rxn\n")

        for sub in substrates:
            i, generated = get_index(sub, metabolites)
            i += MET_START
            if generated:
                node_file.write(f'{i:04d},"{met_file[sub]}",met\n')
            write_met(i, index, edge_file, bool(reversible))

        for prod in products:
            i, generated = get_index(prod, metabolites)
            i += MET_START
            if generated:
                node_file.write(f'{i:04d},"{met_file[prod]}",met\n')
            write_met(index, i, edge_file, bool(reversible))


if __name__ == "__main__":
    total_rxns = []
    cap = len(INPUT_TEXT)
    response_file.write("{\n")
    has_written_first = False
    for i, line in enumerate(INPUT_TEXT[:cap]):
        gene = line[0]
        if gene not in NO_RESPONSE_GENES and gene != "":
            if gene not in cache_file.keys():
                response = requests.get(
                    f"http://bigg.ucsd.edu/api/v2/models/iEK1008/genes/{gene}"
                )
                if response.ok:
                    if not has_written_first:
                        response_file.write(f'  "{gene}": {response.text}\n')
                        has_written_first = True
                    else:
                        response_file.write(f', "{gene}": {response.text}\n')
                    response_rxns = json.loads(response.text)["reactions"]
                else:
                    continue
            else:
                response_rxns = cache_file[gene]["reactions"]

            append_reaction_data(total_rxns, response_rxns)

        print_progress(i + 1, len(INPUT_TEXT[:cap]))
    print()
    response_file.write("}")
    log(f"All reactions: {len(total_rxns)}")
    reaction_file = open(f"{job_name}.rxns", "w")
    reaction_file.write(str(total_rxns))
    reaction_file.close()

    node_file = open(f"{job_name}_node.gdf", "w")
    edge_file = open(f"{job_name}_edge.gdf", "w")
    node_file.write("nodedef>name INTEGER,label VARCHAR,class VARCHAR\n")
    edge_file.write(
        "edgedef>source INTEGER,target INTEGER,weight DOUBLE,directed BOOLEAN\n"
    )

    met_file = json.load(open("biggs_metabolites.json"))
    metabolites = set()
    for i, rxn in enumerate(total_rxns):
        write_reaction(i, rxn, metabolites, node_file, edge_file, met_file)
        print_progress(i + 1, len(total_rxns), 0.01)
    print()
    edge_file.close()
    with open(f"{job_name}_edge.gdf") as edge_file:
        for line in edge_file:
            node_file.write(line)
