# squish tsv files together into one

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), f'..{os.sep}..', 'data') + os.sep

FILE_NAMES = ["Behavioral domains and neural structures.txt",
    "Behavioral domains and their descriptions in the literature.txt",
    "Behavioral Domains and their neural networks.txt",
    "Behavioral Domains, related and  task-related neural regions.txt"]

with open(DATA_DIR+'functional-hirarchy.tsv', 'w') as out:
    
    out.write("### note that after each line of ---'s, there is a row of headers\n")

    for filename in FILE_NAMES:
        title = filename[:-4] # strips .txt
        print(title)
        out.write(f'[{title}]'.center(100,"-")+"\n")
        with open(DATA_DIR+filename, 'r') as f:
            for line in f:
                out.write(line)
        
