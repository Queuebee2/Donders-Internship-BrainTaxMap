from taxmap import DEFAULT_RECORDS_DIR
from pathlib import Path
from distutils.dir_util import copy_tree


def save_list(list, path, headers=[]):
    """Store items as newline separated values in a textfile"""
    with open(path, 'w+') as fout:
        if headers:
            fout.write(",".join(headers) + '\n')
        for item in list:
            fout.write(f"{item}\n")


def load_list(path):
    """Load the single column of values from list of newline separated values in a textfile"""
    items = []
    with open(path) as f:
        for line in f:
            items.append(line.strip('\n'))
    return items


def backup_search(dirname=None):
    from_dir = DEFAULT_RECORDS_DIR
    if not dirname:
        keywords = load_list(Path(DEFAULT_RECORDS_DIR / 'keywords_used.txt'))
        ids = load_list(Path(DEFAULT_RECORDS_DIR / 'ids.txt'))
        dirname = f"{len(keywords)}_kw_{len(ids)}_ids_backup"

    new_dir = Path(f'data/backups/records/{dirname}').resolve()
    new_dir.mkdir(parents=True, exist_ok=True)

    copy_tree(str(from_dir), str(new_dir))
