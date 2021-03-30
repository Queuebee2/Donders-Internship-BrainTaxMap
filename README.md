# Current branch: `cleanup-internship`
Moving on and rework!  

headers with 'goal' are not yet (fully) implemented.

## Different modules to do different things

### `__init__`
Not sure how to use properly
- uses `python-dotenv` to load `DEV_EMAIL` and `NCBI_API_KEY` from `.env`

((almost)independently, hopefully)
### Analysis
Analyse articles
- article reader
- pool of analysis workers
- collect results
- write results to `data/output/...`

### Biosearch
Take lists of keywords and filters to search pubmed and retrieve articles.
Current demo:
```python
from taxmap import biosearch
keywords = ['list', 'of', 'keywords']
mesh_good = ['I', "really", "want", "all", "these", "mesh", "terms"]
mesh_bad = ["no thanks", "to these", "mesh terms"]

unique_ids = biosearch.find_ids(keywords, mesh_good, mesh_bad)
```
`biosearch.find_ids` makes use of `biosearch._buildquery()`

#### finding relevant articles
- create filter
- determine length of query, make sure its not too long including filter. Maybe make number of keywords per batch smaller based on characters used for filters.
- send batches of queries + filters to ncbi `Entrez.esearch`, retrieve article ids
- collect only unique article ids
- store ids under `/data/current/ids.txt`

#### download
in batches of `N` (default 5000), send the found article ids via `Entrez.efetch` and download articles (abstracts, medline format) and write to file (.gzip(!?))
- store articles under `/data/current/abstracts.gzip`


### Visualisers
take output data from `data/output/..` and create graphs
- write graphs to `data/output/graphs`


# BrainTaxMap
The project goal is to map the behavioural taxonomy of the brain by systematically searching and analysing biomedical texts.

# Usage
(in the making)

## goal install
`pip install taxmap`  

## setup 
- Copy and rename `.env.example` -> `.env` and set corresponding values. Most important is `DEV_EMAIL`. `NCBI_API_KEY` is useful to speed up downloads.

## goal demo

### main demo
```python
import taxmap

keywords = taxmap.tools.from_file(keywords_path)
mesh_terms = ['Rodent']
# get ids for lists
ids = taxmap.biosearch.find_ids(keywords, mesh_terms)
# download all relevant articles
taxmap.biosearch.download_abstracts(ids)

# setup and run analyser
analyser = taxmap.analysis.MedlineAnalyser(custom_params)
analyser.run()

# < trim results, process more>

# create graphs
taxmap.visualiser.treemaps...() 

```

### saving/backing-up search results
Running multiple searches? Results from every search can be easily backed up for later. This will contain all files in `data/output/records`.
You can provide a custom directory name for recognition, otherwise, the directory name will be formatted like this: `f"{len(keywords_used)}_kw_{len(ids)}_ids_backup"` and can be found in `data/backups/records/`
```python
from taxmap.tools import backup_search
backup_search()
```
# For the future

# To discuss
- [ ] what should local database be whilst retrieving articles?
    - textfile (large)
    - **(current)** zipped textfile (smaller, max ~40GB for abstracts I think, when downloading whole PubMed which should be avoided)
    - directly into a database ('harder' to set up, makes package less usable)
- [ ] format data for creating treemaps
    - should not contain empty cells -> all empty cells could be 'unidentified'
    - example of skipping cells containing certain string: 
        ```python
        # pandas dataframe
        df = df[~df.stack().str.contains('unidentified').any(level=0)]
        ```
    - todo: how to apply for multiple strings? see https://plotly.com/python/filter/
- [ ] defaults: How do we go about choosing what needs to have a default setting? having many defaults might make tool use easier, also make it harder to adapt tool
- [ ] move these discussions to actual GitHub? 
- [ ] really use `term[Title/Abstract]` to only look in abstracts?
- [ ] `Entrez.read(Entrez.esearch())["idList"]` is a list of [`Bio.Entrez.Parser.StringElement`](https://github.com/biopython/biopython/blob/9e65eebcea1b4d0a497a8f7ebf51c3c1d53cd6e3/Bio/Entrez/Parser.py#L118)
- [ ] nice naming and path for backing up