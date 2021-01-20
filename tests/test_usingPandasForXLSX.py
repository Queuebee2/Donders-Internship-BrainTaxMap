import pandas
import os
from collections import defaultdict
items_of_interest = [
    "Mental, behavioural or neurodevelopmental disorders",
    "Sleep-wake disorders"
]

path = r"\data\lists-other\simpleTabulation.xlsx"

df = pandas.read_excel(path)

print(type(df))

read_state = 0

items_found = defaultdict(set)
current = None

for i, l in enumerate(df['Title']):
    line = l.strip()

    if l.startswith('-') or l.startswith(' ') and read_state:
        # strip nonalphanumeric characters from start of line
        if current is not None:
            while len(line) > 0 and not line[0].isalnum():
                line = line[1:]
            items_found[current].add(line)
    else:
        if line in items_of_interest:
            print(f'At row {i}, found {line} ')
            read_state=1
            current = line
            
        else:
            read_state=0
            current = None
            previous=None


for item in items_of_interest:
    print(item, len(items_found[item]))