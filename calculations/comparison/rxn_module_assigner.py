import json

modules = json.load(open("../../modules/all_modules.json"))
reactions = json.load(open("../../modules/all_reactions.json"))


def find_value(search_key, dictionary):
    for key in dictionary.keys():
        if search_key in dictionary[key]:
            return key


while True:
    reactionID = int(input("Enter reaction: "))
    reactionID = f"R{reactionID:05d}"
    print(reactionID)
    try:
        module_id = find_value(reactionID, reactions)
        try:
            module = find_value(module_id, modules)
            print(module)
        except:
            print("Couldn't find module")
    except:
        print("Couldn't find module id")
