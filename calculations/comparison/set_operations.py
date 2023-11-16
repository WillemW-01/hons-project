import csv
import pyperclip

input_names = [
    "aureus",
    "jcvi",
    "genitalium",
    "minesweeper",
    "pneumoniae",
]

sets = {
    "aureus": set(),
    "jcvi": set(),
    "genitalium": set(),
    "minesweeper": set(),
    "pneumoniae": set(),
}

for path in input_names:
    # input_file = csv.reader(open(f"{path}_rxns.csv"), delimiter="|")
    input_file = open(f"{path}_ecs.txt")
    # next(input_file)
    for line in input_file:
        # sets[path].add(line[0])
        sets[path].add(line.strip())

a = sets["genitalium"].intersection(
    sets["pneumoniae"], sets["aureus"], sets["jcvi"]
)

output = "\n".join(a)
pyperclip.copy(output)
print(output)
print("*copied to clipboard*")
