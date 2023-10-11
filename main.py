import re
import sys
import winsound
import os

import requests
import asyncio
import aiohttp
import time

from Set import Set

START_ID = 2000
OUT_FOLDER = "output"
LOG_FOLDER = "logs"
MAX_THREADS = 52
EDGE_HEADER_1 = (
    "edgedef>source INTEGER,target INTEGER,directed BOOLEAN, weight DOUBLE"
)
EDGE_HEADER_2 = "edgedef>source INTEGER,target INTEGER,label VARCHAR,directed BOOLEAN, weight DOUBLE"
NODE_HEADER = "nodedef>name INTEGER,label VARCHAR,class VARCHAR,color VARCHAR"
THREAD_INPUT_STR = "Enzyme list size is %d. Enter size of batches: "
ERR_ARGUMENT_STR = "ERR: not enough arguments.\nUsage: python main.py <input_file> <job_name> [reaction_file]"

batch_returned = 0
global batch_sent
batch_sent = MAX_THREADS

log_file = open(f"{LOG_FOLDER}/log.log")
global reaction_list
reaction_list = []

# fmt: off
cached_ecs = [
  "2.7.7.6", "2.7.7.7", "3.4.11.18", "6.1.1.21", "2.7.6.1", "2.7.11.1", 
  "3.2.2.23", "3.4.11.18", "3.6.4.13", "4.1.2.13", "4.2.1.11", "6.1.1.6",
  "6.1.1.7", "6.1.1.14", "6.1.1.21", "6.3.4.2"
]
# fmt: on


# fetches the EC numbers from a file in which they are described
def get_enzymes(filepath):
    enzymes = Set()
    with open(filepath) as input_file:
        for line in input_file:
            if "EC:" in line:
                matches = re.findall(
                    r"EC:(\d{1,3}\.\d{1,3}\.\d{1,3}\.[\d]{0,3})", line
                )
                if matches is not None and len(matches) != 0:
                    enzymes.add(matches[0])
    return enzymes


async def get_kegg_result(session: aiohttp.ClientSession, url):
    global batch_returned

    start = url.find("ec:") + 3
    if url[start:] in cached_ecs:
        log(f"Reading file: {url}")
        cached_file = open(f"cached_ecs/{url[start:]}")
        text = cached_file.read()
        batch_returned += 1
        print_progress(batch_returned, batch_sent)
        return 200, str(text), url

    log(f"Request sent: {url}")
    try:
        async with session.get(url + "") as response:
            text = await response.read()
            batch_returned += 1
            print_progress(batch_returned, batch_sent)
            log(f"Response: {url}")
            return response.status, str(text), url
    except Exception as e:
        print(f"Ran into exception with url {url}")
        return None, None, None


async def gather_kegg_results(urls):
    async with aiohttp.ClientSession() as session:
        tasks = [get_kegg_result(session, url) for url in urls]

        start = time.time()
        # wait for all requests to come back
        print_progress(0, batch_sent)
        results = await asyncio.gather(*tasks)
        print()
        end = time.time()
        log(f"Took {end-start:.4} seconds for all requests")

        redo_urls = []
        to_remove = []
        for result in results:
            if len(result) != 3:
                print(f"Error: length of result is not 3: {result}")
                exit(1)

            # unpack status code and the text
            status, text, url = result

            if "Transferred entry" in str(text):
                log(
                    f"Transferred entry for {url}, sending it through to next batch"
                )
                start = text.index("NAME        Transferred to ") + 27
                end = text.index("\\n", start)
                if " and " in text[start:end]:
                    transferred_ecs = text[start:end].split(" and ")
                    for ec in transferred_ecs:
                        redo_urls.append(
                            f"https://rest.kegg.jp/get/ec:{ec.strip()}"
                        )
                else:
                    redo_urls.append(
                        f"https://rest.kegg.jp/get/ec:{text[start:end].strip()}"
                    )
                to_remove.append(result)

        # remove the results that were transferred
        for transferred_result in to_remove:
            log(f"Removing entry from results: {transferred_result[2]}")
            results.remove(transferred_result)

    return results, redo_urls


async def retrieve_all_urls(urls, enzyme_list):
    global batch_returned, batch_sent

    print("- Started retrieval stage")
    results = []
    while len(urls) != 0:
        batch_returned = 0
        batch_sent = len(urls)
        new_results, urls = await gather_kegg_results(urls)
        batch_sent = len(urls)
        results += new_results

        # remember to add the redo urls to the list that weren't there before
        if len(urls) != 0:
            ecs_to_add = [url[url.rfind(":") + 1 :] for url in urls]
            log(f"Adding these ecs to enzyme_list: {ecs_to_add}")
            enzyme_list.add_items(ecs_to_add)
            log(f"Added the ecs")
    print("  - Completed retrieval stage")
    for result in results:
        status, text, url = result
        log(f"{url:38}: {status}")

    log("")

    print("- Started extraction stage")
    total_reactions = []
    for i, result in enumerate(results):
        status, text, url = result
        enzyme = url[url.rfind(":") + 1 :]
        if enzyme == None:
            print(f"Found None enzyme! {url}")
        if int(status) != 200:
            log(f"ERR: look into this url: {url} (status {status})")
            reactions = []
        else:
            reactions = get_reaction_ids(text, enzyme)[1]

        if len(reactions) != 0:
            for rxn in reactions:
                log(f"Reaction list: {reactions}", 1)
                if rxn != "":
                    substrates, products = extract_reaction(rxn)
                    # substrates, products, enzyme
                    total_reactions.append((substrates, products, enzyme))

        print_progress(i + 1, len(results))
    print()
    print("  - Completed extraction stage")

    log("")

    return total_reactions


def get_reaction_ids(text: str, ec_num: str) -> tuple[str, list[str]]:
    log(f"> Getting reaction id for {ec_num}")

    if "Deleted entry" in text:
        log("Look at this enzyme manually", 1)
        return "", []

    start = text.find("ALL_REAC    ")
    if start == -1:
        log("Look at this enzyme manually, no reaction", 1)
        return "", []
    else:
        start += 12

    end = text.find("SUBSTRATE", start)
    if end == -1:
        log("Look at this enzyme manually, no reaction end", 1)
        return "", []

    text = text[start:end]
    text = text.replace("(other)", "")
    text = text.replace("(G)", "")
    text = text.replace("> ", "")
    text = text.replace(";", "")
    text = text.replace("\\n", "")
    text = text.replace("             ", " ")
    text = text.strip()
    reactions = text.split(" ")
    reactions = list(filter(lambda x: x != "", reactions))
    reactions = [rxn.strip() for rxn in reactions]

    return ec_num, reactions


# function to retrieve the data with the KEGG rest api
def extract_reaction(reactionID: str) -> tuple[list[str], list[str]]:
    log(f"> Extracting {reactionID}", 1)

    if (reaction_list != [] and reactionID in reaction_list) or (
        reaction_list == []
    ):
        try:
            response = requests.get(f"https://rest.kegg.jp/get/{reactionID}")
            if not response.ok:
                print(
                    f"Response for {reactionID} was not good: {response.status_code}",
                    2,
                )

            text = response.text

            start = text.find("DEFINITION  ") + 12
            end = text.find("\n", start)
            text = text[start:end]

            reaction = text.split(" <=> ")

            substrates = reaction[0].split(" + ")
            products = reaction[1].split(" + ")

            filter_metabolites(substrates, ["H2O", "H+", "CO2"])
            filter_metabolites(products, ["H2O", "H+", "CO2"])
            log(f"Substrates: {substrates}", 2)
            log(f"Products: {products}", 2)
            return substrates, products
        except Exception as e:
            log(f"Something went wrong: error: {e}")
            return [], []
    else:
        log(f"Reaction is not in the list", 2)
        return [], []


# necessary for writing to the file, since the id has to be a number
def find_id(value, uniq_set: Set):
    return uniq_set.index(value)


# need to change the commas to something else
def fmt_met(metabolite: str):
    return metabolite.replace(",", "#")


def remove_coefficients(substrates: list[str], products: list[str]):
    for i in range(len(substrates)):
        if re.match(r"\d{1,} \w", substrates[i]):
            substrates[i] = substrates[i][substrates[i].find(" ") + 1 :]

    for i in range(len(products)):
        if re.match(r"\d{1,} \w", products[i]):
            products[i] = products[i][products[i].find(" ") + 1 :]


def filter_metabolites(metabolites: list[str], filter_list):
    for filter_item in filter_list:
        if filter_item in metabolites:
            log(f"Removing {filter_item} from {metabolites}", 2)
            metabolites.remove(filter_item)


# fmt: off
def perform_substitutions(metabolites: list[str]):
    for i, metabolite in enumerate(metabolites):
        if metabolite == "NAD(P)+": metabolites[i] = "NADP+"
        if metabolite == "NAD(P)H": metabolites[i] = "NADPH"
        if metabolite.startswith("D-"): metabolites[i] = metabolite[2:]
        if metabolite.startswith("S-"): metabolites[i] = metabolite[2:]
        if metabolite.startswith("(1) "): metabolites[i] = metabolite[4:]
        if metabolite.startswith("(d)"): metabolites[i] = metabolite[3:]
        if metabolite.startswith("(S)-"): metabolites[i] = metabolite[4:]
        if metabolite.startswith("a "): metabolites[i] = metabolite[2:]
        if metabolite.startswith("an "): metabolites[i] = metabolite[3:]
        if metabolite.startswith("beta-"): metabolites[i] = metabolite[5:]
        if metabolite.startswith("alpha-"): metabolites[i] = metabolite[6:]
        if metabolite.startswith("beta-D-"): metabolites[i] = metabolite[7:]
        if metabolite.startswith("alpha-D-"): metabolites[i] = metabolite[8:]
        if "[side 1]" in metabolite: metabolites[i] = metabolite.replace("[side 1]", "")
        if "[side 2]" in metabolite: metabolites[i] = metabolite.replace("[side 2]", "")
        if " (overall reaction)" in metabolite: metabolites[i] = metabolite.replace(" (overall reaction)", "")
        
        metabolites[i] = metabolites[i].strip()
# fmt: on


def update_ec_annotation(enzyme_list):
    # fmt: off
    updated = {
        "2.1.1.":	"2.1.1.185", "2.3.1.":	"2.3.1.12", "1.2.4.":	"2.3.1.61", "2.4.99.":	"2.5.1.145", "2.7.1.":	"2.7.11.1",
        "2.7.7.":	"2.7.7.101", "2.7.8.":	"2.7.8.41", "3.1.21.":	"3.1.21.3", "3.6.5.":	"3.1.3.100", "3.1.3.":	"3.1.3.16",
        "3.4.24.":	"3.4.24.3", "3.6.4.":	"3.6.4.13", "5.4.2.":	"5.4.2.2", "5.4.99.":	"5.4.99.28", "5.99.1.":	"5.6.2.2", "6.1.1.":	"6.1.1.12",
    }  
    # fmt: on

    for i in range(len(enzyme_list)):
        if enzyme_list[i] in updated.keys():
            enzyme_list[i] = updated[enzyme_list[i]]


# function to assist in writing the graph to a file so that it can be passed to gephi
def write_reaction(
    nodes_file,
    edges_file,
    substrates,
    products,
    ec,
    ec_all,
    metabolites,
    has_written,
):
    ec_id = find_id(ec, ec_all)
    if ec_id == None:
        log(f"Couldn't find enzyme {ec}... weird")
    elif not has_written[ec_id]:
        log(f"Writing enzyme {ec_id}")
        nodes_file.write(f"{ec_id},{ec},enzyme\n")
        has_written[ec_id] = True

    remove_coefficients(substrates, products)
    perform_substitutions(substrates)
    perform_substitutions(products)

    i = 0
    for substrate in substrates:
        met_id = metabolites.index(substrate)
        if met_id is None:  # write in nodes file
            met_id = len(metabolites) + START_ID + i
            log(
                f"Did not find metabolite yet, creating id {met_id} (len: {len(metabolites)})"
            )
            nodes_file.write(f"{met_id},{fmt_met(substrate)},metabolite\n")
            i += 1
        else:  # don't write in nodes file
            met_id += START_ID
            log(f"Found metabolite already, has id {met_id}")
            pass

        # write in edges_file
        edges_file.write(f"{met_id},{ec_id},true,1.0\n")

    for product in products:
        met_id = metabolites.index(product)
        if met_id is None:  # write in nodes file
            met_id = len(metabolites) + START_ID + i
            log(
                f"Did not find metabolite yet, creating id {met_id} (len: {len(metabolites)})"
            )
            nodes_file.write(f"{met_id},{fmt_met(product)},metabolite\n")
            i += 1
        else:  # don't write in nodes file
            met_id += START_ID
            log(f"Found metabolite already, has id {met_id}")
            pass

        # write in edges_file
        edges_file.write(f"{ec_id},{met_id},true,1.0\n")

    metabolites.add_items(substrates)
    metabolites.add_items(products)


# prints to stdout how far the program has progressed
def print_progress(i, length):
    completed = i / length
    not_completed = 80 - int(80 * completed)
    # print(f"{completed:.2f}% complete")
    print(
        f"[{'#' * int(80 * completed)}{' ' * (not_completed)}] {completed * 100:.2f}",
        end="\r",
    )
    time.sleep(0.01)
    pass


# writes to the log file what has happened
def log(message: str, level=0):
    global log_file
    indent = "  " * level
    log_file.write(f"{indent}{message}\n")
    log_file.flush()
    os.fsync(log_file.fileno())


if __name__ == "__main__":
    start = time.time()

    if len(sys.argv) < 3:
        exit(ERR_ARGUMENT_STR)

    input_file = sys.argv[1]
    job_name = sys.argv[2]
    if len(sys.argv) == 4:
        reaction_input = open(sys.argv[3])
        reaction_list = reaction_input.readlines()
        reaction_list = [rxn.strip() for rxn in reaction_list]
        reaction_input.close()
    else:
        reaction_list = []

    log_file = open(f"{LOG_FOLDER}/{job_name}.log", "w")
    nodes_file = open(f"{OUT_FOLDER}/{sys.argv[2]}_nodes.gdf", "w")
    edges_file = open(f"{OUT_FOLDER}/{sys.argv[2]}_edges.gdf", "w")

    # fmt: off
    if sys.argv[1] != "default":
        enzyme_list = get_enzymes(sys.argv[1])
    else:
        enzyme_list = Set(items=["6.1.1.21","2.7.4.9","2.5.1.19","1.2.1.11","2.5.1.1"])
    # fmt: on

    log(f"Before updating ecs: {str(enzyme_list)}")
    update_ec_annotation(enzyme_list)
    log(f"After updating ecs: {str(enzyme_list)}")

    thread_input = input(THREAD_INPUT_STR % len(enzyme_list))
    if thread_input != "":
        MAX_THREADS = int(thread_input)

    nodes_file.write(f"{NODE_HEADER}\n")
    edges_file.write(f"{EDGE_HEADER_1}\n")

    urls = [f"https://rest.kegg.jp/get/ec:{ec_num}" for ec_num in enzyme_list]

    if len(urls) > MAX_THREADS:
        reactions = []
        k = len(urls) / MAX_THREADS
        if k - int(k) > 0:  # if there is a fraction
            k = int(k) + 1
        else:
            k = int(k)

        for i in range(int(k)):
            print(f"Sent out batch {i+1} / {k}")
            urls_batch = urls[MAX_THREADS * i : MAX_THREADS * (i + 1)]
            reactions += asyncio.run(retrieve_all_urls(urls_batch, enzyme_list))
            batch_returned = 0
            print(f"Got back batch {i + 1} / {k}")

            # prevents being rate limited by sleeping for 2 minutes
            if int(k) > 2 and i < k - 1:
                time.sleep(120)

    else:
        reactions = asyncio.run(retrieve_all_urls(urls, enzyme_list))

    print("- Starting writing stage")

    counter = 0
    metabolites = Set()
    has_written = [False for i in range(len(enzyme_list))]
    for reaction in reactions:
        if all(reaction):
            substrates, products, enzyme = reaction
            log("Writing metabolites to graph")

            write_reaction(
                nodes_file,
                edges_file,
                substrates,
                products,
                enzyme,
                enzyme_list,
                metabolites,
                has_written,
            )

        print_progress(counter, len(reactions))

        log("")
        counter += 1
    print_progress(counter, len(reactions))
    print()

    edges_file.close()

    # appends the edges file to the nodes file
    with open(f"{OUT_FOLDER}/{sys.argv[2]}_edges.gdf") as edges_file:
        for line in edges_file:
            nodes_file.write(line)

    nodes_file.close()
    print("  - Completed writing stage")

    log("<== Finished ==>")
    print(f"Took {time.time() - start:.4f} seconds to run.")
    log(f"Took {time.time() - start:.4f} seconds to run.")
    log_file.close()

    # plays a notification sound when finished
    winsound.Beep(1500, 500)
    winsound.Beep(1500, 200)
    winsound.Beep(1500, 500)
