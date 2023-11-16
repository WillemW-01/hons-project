import asyncio
import json
import os
import re
import sys
import time

# non-stdlib dependencies
import aiohttp
import requests

# custom set object
from Set import Set

# global constants
START_ID = 2000
OUT_FOLDER = "output"
LOG_FOLDER = "logs"
MAX_THREADS = 52
EDGE_HEADER_1 = "edgedef>source INTEGER,target INTEGER,directed BOOLEAN,weight DOUBLE"
EDGE_HEADER_2 = (
    "edgedef>source INTEGER,target INTEGER,label VARCHAR,directed BOOLEAN,weight DOUBLE"
)
NODE_HEADER = (
    "nodedef>name INTEGER,label VARCHAR,class VARCHAR,module VARCHAR,color VARCHAR"
)
THREAD_INPUT_STR = "Enzyme list size is %d. Enter size of batches: "
USAGE_STR = "ERR: not enough arguments.\nUsage: python main.py <input_file> <job_name> [reaction_file]"

# fmt: off
MODULE_COLORS = json.load(open("modules/module_colors.json"))
CACHED_REACTIONS = json.load(open("modules/cached_reactions.json"))
CACHED_REACTION_IDS = json.load(open("modules/cached_reaction_ids.json"))
CACHED_PATHWAYS = json.load(open("modules/map_conversion.json"))
CACHED_ECS = [
  "1.5.1.3", "2.1.1.198", "2.3.1.234", "2.4.2.1", "2.4.2.9", "2.5.1.6", "2.5.1.145", "2.7.1.11",
  "2.7.1.23", "2.7.1.30", "2.7.1.40", "2.7.1.194", "2.7.4.8", "2.7.6.1", "2.7.7.2", "2.7.7.6",
  "2.7.7.7", "2.7.7.18", "2.7.8.5", "2.7.8.7", "2.7.11.1", "3.1.1.29", "3.1.3.7", "3.1.21.10",
  "3.2.2.23", "3.4.11.18", "3.4.21.53", "3.4.23.36", "3.5.1.88", "3.6.4.13", "4.1.2.13", "4.2.1.11",
  "5.1.3.2", "5.3.1.9", "5.4.2.2", "5.6.2.1", "5.6.2.2", "6.1.1.2", "6.1.1.6", "6.1.1.7",
  "6.1.1.10", "6.1.1.11", "6.1.1.12", "6.1.1.14", "6.1.1.20", "6.1.1.21", "6.3.4.2", "6.3.4.21",
  "7.1.1.2"
]
OVERVIEW_MAPS = [
  "ec01100", "ec01110", "ec01120", "ec01200", "ec01210", "ec01212",
  "ec01230", "ec01232", "ec01250", "ec01240", "ec01220",
]
# fmt: on

# global variables
global batch_sent
global reaction_list
enzyme_modules = {}
reaction_list = []
batch_returned = 0
batch_sent = MAX_THREADS
log_file = open(f"{LOG_FOLDER}/log.log")


# fetches the EC numbers from a file in which they are described
def get_enzymes(filepath: str):
    """
    Extract the EC numbers from a file.

    The function searches line by line, and searches for the pattern
    `EC:a.b.c.d` where a,b,c,d are numbers.
    It adds matches of this pattern to a custom set object.
    Only EC numbers with at least 3 numbers are retrieved (e.g. 1.8.4.-)

    :param filepath: the path of the file to read
    :return: a set of non-duplicate EC numbers
    """
    enzymes = Set()
    with open(filepath) as input_file:
        for line in input_file:
            if "EC:" in line:
                matches = re.findall(r"EC:(\d{1,3}\.\d{1,3}\.\d{1,3}\.[\d]{1,3})", line)
                if matches is not None and len(matches) != 0:
                    enzymes.add(matches[0])
    return enzymes


async def get_kegg_result(session: aiohttp.ClientSession, url: str) -> tuple:
    """
    Retrieves a KEGG EC entry text.

    The function first looks if the url is specified in the list of cached
    EC numbers (entries that take very long to retrieve), and will return
    its contents if it is true.

    If not, it will send a GET request to KEGG and return the response text.

    For both cases it will update the progress bar being drawn to the terminal.

    Will return a tuple of Nones if there was an error along the way.

    :param session: the session created by the atiohttp package that is used to run these requests concurrently.
    :param url: the url that is used in the GET request.
    :return: a tuple containing the response code, response text, and the original url.

    """
    global batch_returned

    start = url.find("ec:") + 3  # isolates the EC number from the url
    if url[start:] in CACHED_ECS:
        log(f"Reading file: {url}")
        cached_file = open(f"cached_ecs/{url[start:]}.txt")
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
        print(f"Exception: {e}")
        return None, None, None


async def gather_kegg_results(urls: list[str]) -> tuple[list[tuple], list[str]]:
    """
    Sends requests to the KEGG API concurrently.

    This function will send requests to KEGG with all the urls contained in the
    input list. It will use the aiohttp library to handle the session, and the
    asyncio library to send all the requests together.

    It also calculates the time for the batch of requests.

    After all the responses are returned, it determines whether some of the EC
    numbers were transferred to different entries. It adds these urls to a list
    that will be redone in the next iteration called by `retrieve_all_urls`.

    If it was not transferred, then the results are simply appended to a list.
    This list is returned.

    :param urls: the list of urls to be sent to KEGG
    :return: a tuple of a list of results, and a list of urls that were transferred.
    """
    async with aiohttp.ClientSession() as session:
        tasks = [get_kegg_result(session, url) for url in urls]

        start = time.time()
        print_progress(0, batch_sent)

        # wait for all requests to come back
        results = await asyncio.gather(*tasks)

        print()
        end = time.time()
        log(f"Took {end-start:.4} seconds for all requests")

        redo_urls = []
        to_remove = []
        for result in results:
            if len(result) != 3:  # this should not happen
                print(f"Error: length of result is not 3: {result}")
                exit(1)

            # unpack status code and the text
            status, text, url = result

            # will redo this url if the EC entry was transferred to any other
            # number of entries
            if "Transferred entry" in str(text):
                log(f"Transferred entry for {url}, sending it through to next batch")
                start = text.index("NAME        Transferred to ") + 27
                end = text.index("\\n", start)
                if " and " in text[start:end]:
                    transferred_ecs = text[start:end].split(" and ")
                    for ec in transferred_ecs:
                        redo_urls.append(f"https://rest.kegg.jp/get/ec:{ec.strip()}")
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


def get_reaction_ids(text: str, ec_num: str) -> tuple[str, list[str]]:
    """
    Extracts the reaction identifier(s) from a given KEGG EC entry.

    This will essentially retrieve all the reactions that are catalysed by the
    enzyme specified by ec_num.

    Some characters that were added to reaction identifiers in order to clean
    the list are removed.

    :param text: the response text for the given KEGG EC entry.
    :param ec_num: the EC number of the response given.
    :return: a tuple containing the input EC number, as well as all the reaction
    identifiers.
    """

    log(f"> Getting reaction id for {ec_num}", 1)

    if ec_num in CACHED_REACTION_IDS.keys():
        return ec_num, CACHED_REACTION_IDS[ec_num]

    if "Deleted entry" in text:
        log("Look at this enzyme manually, deleted enzyme", 2)
        return "", []

    start = text.find("ALL_REAC    ")
    if start == -1:
        log("Look at this enzyme manually, no reaction", 2)
        return "", []
    else:
        start += 12

    end = text.find("SUBSTRATE", start)
    if end == -1:
        log("Look at this enzyme manually, no reaction end", 2)
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


def get_pathway_ids(text: str, ec_num: str):
    log(f"> Getting pathway ids for {ec_num}", 1)
    should_add = False
    pathways = []
    split_string = "\\n" if "\\n" in text else "\n"
    for line in text.split(split_string):
        if line.startswith("PATHWAY"):
            should_add = True
        if (
            line.startswith("ORTHOLOGY")
            or line.startswith("GENES")
            or line.startswith("DBLINKS")
        ):
            should_add = False

        if should_add:
            matches = re.findall(r"ec\d{5}", line)
            if matches is not None and matches[0] not in OVERVIEW_MAPS:
                log(f"Adding pathway {matches[0]}", 2)
                pathways.append(matches[0])

    return pathways


def get_module_ids(pathways: list[str]):
    log(f"> Getting module ids", 1)
    module_ids = []
    for pathway in pathways:
        log(f"> Finding modules for {pathway}", 2)

        if pathway in CACHED_PATHWAYS:
            log(f"Modules were cached. Found: {CACHED_PATHWAYS[pathway]}", 3)
            module_ids += CACHED_PATHWAYS[pathway]
            continue

        response = requests.get(f"https://rest.kegg.jp/get/{pathway}")
        if response.ok:
            text = response.text
            should_add = False
            for line in text.split("\n"):
                if line.startswith("MODULE"):
                    should_add = True
                if (
                    line.startswith("DBLINKS")
                    or line.startswith("ENZYME")
                    or line.startswith("REL_PATHWAY")
                    or line.startswith("COMPOUND")
                ):
                    should_add = False

                if should_add:
                    matches = re.findall(r"M\d{5}", line)
                    if matches is not None:
                        log(f"Adding module id {matches[0]}", 3)
                        module_ids.append(matches[0])
    return module_ids


def extract_reaction(reactionID: str) -> tuple[list[str], list[str]]:
    """
    Extracts substrates and products from a KEGG reaction entry.

    This function will send a request to the KEGG REACTION database for the given
    reaction identifier. It first checks if the reaction should be analysed, by
    first checking that the reactionID is contained in the optional reaction list
    specified in program input.

    If it should be analysed, it will separate the reaction definition by the
    reversibility sign, and split the substrates and products by the plus sign.
    It will remove the inorganic metabolites water, proton, and carbon dioxide,
    as these should not be involved in the metabolic representation.

    It then returns the substrates and products.

    :param reactionID: the reaction identifier to be analysed.
    :return a tuple of the substrates and products, which are in themselves lists of strings.

    """
    log(f"> Extracting {reactionID}", 1)

    if (reaction_list != [] and reactionID in reaction_list) or (reaction_list == []):
        if reactionID in CACHED_REACTIONS.keys():
            substrates, products = (
                CACHED_REACTIONS[reactionID]["Substrates"],
                CACHED_REACTIONS[reactionID]["Products"],
            )
            log(f"Substrates: {substrates}", 2)
            log(f"Products: {products}", 2)
            return substrates, products
        try:
            response = requests.get(f"https://rest.kegg.jp/get/{reactionID}")
            if not response.ok:
                log(
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

            # these should not be in metabolic graph, removes them
            # for every x in substrates and products, only retains it if it is
            # not water, a proton or carbon dioxide.
            filter_list = ["H2O", "H+", "CO2", "H20"]
            substrates = list(filter(lambda x: x not in filter_list, substrates))
            products = list(filter(lambda x: x not in filter_list, products))

            log(f"Substrates: {substrates}", 2)
            log(f"Products: {products}", 2)

            return substrates, products
        except Exception as e:
            log(f"Something went wrong: error: {e}")
            return [], []
    else:
        log(f"Reaction is not in the list", 2)
        return [], []


def assign_modules(enzyme, modules):
    global enzyme_modules

    log(f"> Assigning modules to {enzyme}", 1)

    # reaction_file = open("modules/aureus_reactions.json")
    # reactions_dict = json.load(reaction_file)
    # module_id = ""
    # for key in reactions_dict.keys():
    #     if rxn_id in reactions_dict[key]:
    #         module_id = key
    #         break

    if enzyme not in enzyme_modules.keys():
        log(f"Enzyme not in the list yet, creating a new set.", 2)
        enzyme_modules[enzyme] = {}

    if len(modules) == 0:
        return

    modules_file = open("modules/all_modules.json")
    modules_dict = json.load(modules_file)

    # for each of the found module ids, find the corresponding module
    # e.g. M00019 -> Amino acid metabolism
    for module in modules:
        # for each of the possible modules (e.g. Carbohydrate metabolism, etc)
        for key in modules_dict.keys():
            if module in modules_dict[key]:
                if key not in enzyme_modules[enzyme].keys():
                    enzyme_modules[enzyme][key] = 1
                else:
                    enzyme_modules[enzyme][key] += 1

    if len(enzyme_modules[enzyme]) != 0:
        log(f"Modules assigned to {enzyme}: {enzyme_modules[enzyme]}", 2)
    else:
        log(f"No modules assigned to {enzyme}")


def assign_reactions(enzyme, reaction_ids):
    global enzyme_modules

    log(f"> Assigning modules to {enzyme}", 1)

    # reaction_file = open("modules/aureus_reactions.json")
    # reactions_dict = json.load(reaction_file)
    # module_id = ""
    # for key in reactions_dict.keys():
    #     if rxn_id in reactions_dict[key]:
    #         module_id = key
    #         break

    if enzyme not in enzyme_modules.keys():
        log(f"Enzyme not in the list yet, creating a new set.", 2)
        enzyme_modules[enzyme] = {}

    if len(reaction_ids) == 0:
        return

    modules_file = open("modules/all_reactions_new.json")
    modules_dict = json.load(modules_file)

    # for each of the found module ids, find the corresponding module
    # e.g. M00019 -> Amino acid metabolism
    for rxn in reaction_ids:
        # for each of the possible modules (e.g. Carbohydrate metabolism, etc)
        for key in modules_dict.keys():
            if rxn in modules_dict[key]:
                if key not in enzyme_modules[enzyme].keys():
                    enzyme_modules[enzyme][key] = 1
                else:
                    enzyme_modules[enzyme][key] += 1

    if len(enzyme_modules[enzyme]) != 0:
        log(f"Modules assigned to {enzyme}: {enzyme_modules[enzyme]}", 2)
    else:
        log(f"No modules assigned to {enzyme}")


async def process_all_urls(urls: list[str], enzyme_list: Set) -> list[tuple]:
    """
    Retrieves all reactions catalysed by the list of EC numbers in the urls.

    Sends all the urls to KEGG and retrieves the EC entries. If EC entries were
    transferred, it will redo the transferred url. For each EC entry, it will
    retrieve the corresponding reaction identifier(s). For each reaction,
    it will extract the substrates and products and add these to a list. After
    every reaction was processed, the list of reactions are returned.

    :param urls: the list of urls to be sent to KEGG
    :enzyme_list: the global list of EC numbers that must be updated for transferred entries.
    :return: a list containing tuples of (substrates, products, ecs) for each reaction.
    """
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

    print("  - Completed retrieval stage")

    for result in results:
        status, text, url = result
        log(f"{url:38}: {status}")
    log("")

    print("- Started extraction stage")
    total_reactions = []
    for i, result in enumerate(results):
        status, text, url = result
        enzyme = url[url.rfind(":") + 1 :]  # extracts EC from url

        if enzyme == None:
            print(f"Found None enzyme! {url}")  # this should not happen
        if int(status) != 200:  # http error along the way, should discard
            log(f"ERR: look into this url: {url} (status {status})")
            reactions = []
        else:
            log(f"> Working with {enzyme}")
            # pathways = get_pathway_ids(text, enzyme)
            # modules = get_module_ids(pathways)
            # assign_modules(enzyme, modules)

            reactions = get_reaction_ids(text, enzyme)[1]
            assign_reactions(enzyme, reactions)

        if len(reactions) != 0:  # if the list of reactions contained something
            log(f"Reaction list: {reactions}", 1)
            for rxn in reactions:
                if rxn != "":
                    substrates, products = extract_reaction(rxn)

                    total_reactions.append((substrates, products, enzyme))

        print_progress(i + 1, len(results))

    print()
    print("  - Completed extraction stage")

    log("")

    return total_reactions


def fmt_met(metabolite: str):
    """This function simply replaces commas with a hastag in order to make the
    output file a csv. For example: 1,6-fructose bisphosphate will break the csv format.

    :param metabolite: the metabolite to format.
    :return the input string with the comma replaced.
    """
    return metabolite.replace(",", "#")


def remove_coefficients(substrates: list[str], products: list[str]):
    """
    This function removes the coefficients before metabolites in the metabolic
    reaction. For example `(2) folate` becomes `folate`. It does this for both
    substrates and products.

    The function works by matching the pattern `(d) w` where d is a number and is
    contained within parentheses and w is the metabolite after the coefficient.
    The function will retain everything after the space between the coefficient
    and the metabolite if such a match was found.

    :param substrates: the list of substrate strings in the current reaction.
    :param products: the list of product strings in the current reaction.
    :return None, however the lists were implictly updated.
    """

    for i in range(len(substrates)):
        if re.match(r"\d{1,} \w", substrates[i]):
            substrates[i] = substrates[i][substrates[i].find(" ") + 1 :]

    for i in range(len(products)):
        if re.match(r"\d{1,} \w", products[i]):
            products[i] = products[i][products[i].find(" ") + 1 :]


# fmt: off
def perform_substitutions(metabolites: list[str]):
    """
    This function makes the metabolite labels easier to interprate.
    
    This is done by ensuring that there are no "pseudo"-duplicate entries,
    where two metabolites don't share the same label, but are in fact the same
    metabolite. For example "alpha-D-glucose" and "D-glucose"
    (alpha configuration assumed to be default).
    
    Also removes some text in KEGG annotations that are unnecessary for further
    analysis.
    
    :param metabolites: the list of metabolites to be changed.
    :return: None, however the metabolite list was updated implicitly. 
    """

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


def update_ec_annotation(enzyme_list: Set):
    """
    This function will update partial EC numbers (e.g. 1.8.4.-) to an appropriate
    EC number that was determined manually.

    :param enzyme_list: the global set of EC numbers that is being analysed.
    :return: None, however, the input Set object was changed implicitly.
    """

    # fmt: off
    updated = {
        "2.1.1.":	"2.1.1.185", "2.3.1.":	"2.3.1.12", "1.2.4.":	"2.3.1.61",
        "2.4.99.":	"2.5.1.145", "2.7.1.":	"2.7.11.1", "2.7.7.":	"2.7.7.101",
        "2.7.8.":	"2.7.8.41", "3.1.21.":	"3.1.21.3", "3.6.5.":	"3.1.3.100",
        "3.1.3.":	"3.1.3.16", "3.4.24.":	"3.4.24.3", "3.6.4.":	"3.6.4.13",
        "5.4.2.":	"5.4.2.2", "5.4.99.":	"5.4.99.28", "5.99.1.":	"5.6.2.2",
        "6.1.1.":	"6.1.1.12",
    }  
    # fmt: on

    for i in range(len(enzyme_list)):
        if enzyme_list[i] in updated.keys():
            enzyme_list[i] = updated[enzyme_list[i]]


def get_one_module(ec):
    module_list = enzyme_modules[ec]
    # priority_list = [
    #     "Carbohydrate metabolism",
    #     "Energy metabolism",
    #     "Lipid metabolism",
    #     "Nucleotide metabolism",
    #     "Amino acid metabolism",
    #     "Glycan metabolism",
    #     "Metabolism of cofactors and vitamins",
    #     "Biosynthesis of terpenoids and polyketides",
    #     "Biosynthesis of other secondary metabolites",
    #     "Xenobiotics biodegradation",
    # ]
    if len(module_list) == 0:
        return None
    if len(module_list) == 1:
        return list(module_list.keys())[0]
    elif max(module_list.values()) == 1:
        priority_list = [
            "Carbohydrate metabolism",
            "Energy metabolism",
            "Lipid metabolism",
            "Nucleotide metabolism",
            "Amino acid metabolism",
            "Metabolism of cofactors and vitamins",
            "Biosynthesis of terpenoids and polyketides",
            "Biosynthesis of other secondary metabolites",
            "Pathogenicity",
        ]

        for module in priority_list:
            if module in module_list.keys():
                return module
    else:
        max_count = max(module_list.values())
        for key in module_list.keys():
            if module_list[key] == max_count:
                return key

    return None


def write_reaction(
    nodes_file,
    edges_file,
    substrates: list[str],
    products: list[str],
    ec: str,
    ec_all: Set,
    metabolites: Set,
    has_written: list[bool],
):
    """
    This function writes the reaction given by the tuple (substrates, products, ec) to
    the output file.

    It does this by first writing the EC number as a node if it has not been written yet.
    Then it cleans up the labels of the metabolites. It will then fetch an ID
    for each metabolite, and create a new one if the metabolite has not been
    written to the file yet. These IDs are important for Gephi (visualisation program)
    to function.

    Then, for each substate, it will write an edge between the substrate and the
    EC catalysing the reaction. Similarly, for each product, an edge between the
    EC and the product will be written.

    The substrates and products of this reaction is then added to a list that is
    used to determine their ID for further calls of this function.

    :param nodes_file: the file where the nodes of the graph are written
    :param edges_file: the file where the edges of the graph are written
    :param substrates: the list of substrates in the current reaction
    :param products: the list of products in the current reaction
    :param ec: the EC of the enzyme catalysing the current reaction
    :param ec_all: the global list of enzyme numbers, used to determine the id of the current EC
    :param metabolites: the list of all metabolites seen so far, used to determine the ids of substrates and products
    :param has_written: used to keep track if a certain EC number was already written as a node.
    :return: None
    """
    # ec_id = find_id(ec, ec_all)
    ec_id = ec_all.index(ec)
    if ec_id == None:
        log(f"Couldn't find enzyme {ec}... weird")
    elif not has_written[ec_id]:
        log(f"Writing enzyme {ec_id}")
        module = get_one_module(ec)
        if module is None:
            nodes_file.write(f"{ec_id},{ec},enzyme,,{MODULE_COLORS['None']}\n")
        else:
            nodes_file.write(f"{ec_id},{ec},enzyme,{module},{MODULE_COLORS[module]}\n")
        has_written[ec_id] = True

    remove_coefficients(substrates, products)
    perform_substitutions(substrates)
    perform_substitutions(products)

    i = 0
    for substrate in substrates:
        met_id = metabolites.index(substrate)
        if met_id is None:  # write in nodes file, haven't found yet
            met_id = len(metabolites) + START_ID + i
            log(
                f"Did not find metabolite yet, creating id {met_id} (len: {len(metabolites)})"
            )
            nodes_file.write(
                f"{met_id},{fmt_met(substrate)},metabolite,,{MODULE_COLORS['None']}\n"
            )
            i += 1
        else:  # don't write in nodes file, already has id
            met_id += START_ID
            log(f"Found metabolite already, has id {met_id}")
            pass

        # write in edges_file
        edges_file.write(f"{met_id},{ec_id},true,1.0\n")

    metabolites.add_items(substrates)

    i = 0
    for product in products:
        met_id = metabolites.index(product)
        if met_id is None:  # write in nodes file, haven't found yet
            met_id = len(metabolites) + START_ID + i
            log(
                f"Did not find metabolite yet, creating id {met_id} (len: {len(metabolites)})"
            )
            nodes_file.write(
                f"{met_id},{fmt_met(product)},metabolite,,{MODULE_COLORS['None']}\n"
            )
            i += 1
        else:  # don't write in nodes file, found already
            met_id += START_ID
            log(f"Found metabolite already, has id {met_id}")
            pass

        # write in edges_file
        edges_file.write(f"{ec_id},{met_id},true,1.0\n")

    metabolites.add_items(products)


def print_progress(i, length):
    """
    Updates the progress bar being drawn in the terminal.

    This function simply draws a line of "#", scaled to 80, depending on how
    far the progress has been (i / length). It also prints out a carriage return
    such that the progress bar overwrites itself.

    It also sleeps for 10 msec.

    :param i: the progress made so far (assumed 0 < i <= length)
    :param length: the total amount of progress that can be made
    :return: None
    """

    completed = i / length
    not_completed = 80 - int(80 * completed)
    hash_part = "#" * int(80 * completed)
    empty_part = " " * (not_completed)

    print(f"[{hash_part}{empty_part}] {completed * 100:.2f}%", end="\r")
    time.sleep(0.01)


def log(message: str, level=0):
    """
    Writes a string to the log file and not to stdout.

    The function also includes a level, which will simply append two spaces
    before the input string. The log_file is a global TextIOWrapper object that
    is created during runtime.

    The log file is also flushed so that the output can be visualized in real time.

    :param message: the string to write to the log file
    :param level: the amount of indendation for the string
    :return: None
    """

    # global log_file
    indent = "  " * level
    log_file.write(f"{indent}{message}\n")
    log_file.flush()
    os.fsync(log_file.fileno())


def report_time(start):
    end = time.time()
    duration = end - start
    if duration > 120:
        minutes = int(duration // 60)
        seconds = int(duration - minutes * 60)
        print(f"Took {minutes:d} minutes and {seconds:d} seconds to run.")
        log(f"Took {minutes:d} and {seconds:d} seconds to run.")
    else:
        print(f"Took {duration:.4f} seconds to run.")
        log(f"Took {duration:.4f} seconds to run.")


def main():
    """
    The main function to execute the program based on the input parameters. Will
    produce the graph file based on the reactions specified by the input EC numbers.
    """
    global MAX_THREADS, log_file, enzyme_modules, reaction_list

    start = time.time()  # start tracking the execution time

    if len(sys.argv) < 3:  # did not run correctly, print out usage string
        exit(USAGE_STR)

    ec_input_file = sys.argv[1]
    job_name = sys.argv[2]

    if len(sys.argv) == 4:  # if a reaction file was specified
        reaction_input = open(sys.argv[3])
        reaction_list = reaction_input.readlines()
        reaction_list = [rxn.strip() for rxn in reaction_list]
        reaction_input.close()
    else:
        reaction_list = []

    # opens the files used for output
    log_file = open(f"{LOG_FOLDER}/{job_name}.log", "w")
    nodes_file = open(f"{OUT_FOLDER}/{sys.argv[2]}_nodes.gdf", "w")
    edges_file = open(f"{OUT_FOLDER}/{sys.argv[2]}_edges.gdf", "w")

    # fmt: off
    # if no EC input file was specified, uses a testing set of ECs.
    # otherwise, fetches the ECs from the input file
    if sys.argv[1] != "default":
        enzyme_list = get_enzymes(ec_input_file)
    else:
        enzyme_list = Set(items=["6.1.1.21","2.7.4.9","2.5.1.19","1.2.1.11","2.5.1.1"])
    # fmt: on

    # Updates partial EC annotations
    log(f"Before updating ecs: {str(enzyme_list)}")
    update_ec_annotation(enzyme_list)
    log(f"After updating ecs: {str(enzyme_list)}")

    # gets the size of the batches to be sent to KEGG.
    # A low enough value prevents being rate limited (~50)
    thread_input = input(THREAD_INPUT_STR % len(enzyme_list))
    if thread_input != "":
        MAX_THREADS = int(thread_input)

    urls = [f"https://rest.kegg.jp/get/ec:{ec_num}" for ec_num in enzyme_list]

    if len(urls) > MAX_THREADS:  # if the ECs must be separated into batches
        reactions = []
        k = len(urls) / MAX_THREADS  # k = number of batches
        if k - int(k) > 0:  # if there is a fraction
            k = int(k) + 1  # e.g. 3.5 -> 4
        else:
            k = int(k)

        for i in range(int(k)):  # for every batch i
            print(f"Sent out batch {i+1} / {k}")
            urls_batch = urls[MAX_THREADS * i : MAX_THREADS * (i + 1)]
            reactions += asyncio.run(process_all_urls(urls_batch, enzyme_list))
            batch_returned = 0  # update global variable for print_progress
            print(f"Got back batch {i + 1} / {k}")

            # prevents being rate limited by sleeping for 2 minutes
            if int(k) > 2 and i < k - 1:
                time.sleep(120)
    else:  # if the number of ECs is small enough to not be broken up
        reactions = asyncio.run(process_all_urls(urls, enzyme_list))

    print("- Starting writing stage")

    nodes_file.write(f"{NODE_HEADER}\n")
    edges_file.write(f"{EDGE_HEADER_1}\n")

    counter = 0
    metabolites = Set()
    has_written = [False for i in range(len(enzyme_list))]
    for reaction in reactions:
        if all(reaction):
            substrates, products, enzyme = reaction
            log("Writing metabolites to graph")

            # fmt: off
            write_reaction(nodes_file, edges_file, substrates, products, enzyme,
                enzyme_list, metabolites, has_written)  # fmt: on

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
    os.remove(f"{OUT_FOLDER}/{sys.argv[2]}_edges.gdf")  # to clean up output folder
    print("  - Completed writing stage")

    log("<== Finished ==>")
    report_time(start)
    log_file.close()


if __name__ == "__main__":
    main()
