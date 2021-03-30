import logging
import sys
from pathlib import Path
from taxmap import tools


logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(name)s -- %(levelname)s : %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)


def demo():
    from taxmap import biosearch
    from taxmap import tools
    INPUT_PATH = Path('data/input').resolve()
    OUTPUT_PATH = Path('data/output').resolve()

    # to load keywords from file
    # search_keywords = tools.load_list(Path(INPUT_PATH/'pruned_structures_list.txt'))
    
    # manually set keywords
    search_keywords =[ "Striatum", "Somatosensory cortex", "visual Cortex", "Amygdala", "cerebellum", "prefrontal cortex", "premotor cortex"]

    mesh_include = ["Mice"]
    mesh_exclude = ["Human"]

    ids = biosearch.find_ids(
        search_keywords, 
        mesh_include=mesh_include, 
        mesh_exclude=mesh_exclude)
    
    ids = [int(x) for x in ids]

    
    # comment out top method if want to use saved ids
    # ids = tools.load_list(Path(OUTPUT_PATH/'records'/'ids.txt'))
    print(f"FOUND {len(ids)} ids!")
    biosearch.download_abstracts(ids)


if __name__ == '__main__':
    print('running the demo')
    demo()