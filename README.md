Donders-Internship-BrainTaxMap

# Project description  
Map brain structure to function based on literature  
  
# Current features  
- parse JSON structural hirachry (`root:[{child:[{child}]]`)  
- parse dotnotation `Cognition.Language.Orthography` -> 3 levels of 'functional hierarchy'  
- search pubmed  
    - currently only select [PMC](https://www.ncbi.nlm.nih.gov/pmc/) articles  
    - amount of results can be set  
    - searched keywords are stored to prevent re-searching, can be reset 
    - one issue where database connection closes, dirty solution used at this moment
- store data in [Neo4j](https://neo4j.com)  
    - visualisations in neo4j bloom and access to 'playground'  
    - current stats:  
        - 45k nodes 
        - 75k relationships  


# Notes on `Neo4j`, `graphs`, `py2neo`, `gephi`
work in progress

### NEO4J Browser commands
- Find all nodes of one type:   
 `MATCH (n:brainstructure) RETURN n`
- Find ALL nodes  
 `Match (n)-[r]->(m) Return n,r,m`
- See Scheme  
`:scheme`
- See constraints  
`CALL db.constraints`
- Find relations between specific nodes  
`Match (n:nodelabel)--(m:othernodelabel) return n,m`
- Find specific relations  
`MATCH (n1)-[r:SPECIFIC_RELATION]->(n2) RETURN r, n1, n2`
- Replace relations between specific nodes  
```
MATCH (n1:label)-[old:OLD_RELATION]->(n2:label) 
MERGE (n1)-[new:NEW_RELATION]->(n2)
DELETE old
```
- [more here](https://www.remwebdevelopment.com/blog/sql/some-basic-and-useful-cypher-queries-for-neo4j-201.html)



### RESOURCES  
- [batched updates](https://medium.com/neo4j/5-tips-tricks-for-fast-batched-updates-of-graph-structures-with-neo4j-and-cypher-73c7f693c8cc)
- [hosting  through azure](https://azuremarketplace.microsoft.com/en-us/marketplace/apps?search=neo4j&page=1)
- [Medline Records Object Attributes](https://biopython.org/docs/1.75/api/Bio.Medline.html)
    
#### TODO  
- bhttps://www.biostars.org/p/89478/ extract citations from PMC articles for art->[CITED]->art relationship  
- look in to Chord dependency diagram or Hierarchical Edge Bundling or [force layout](https://github.com/d3/d3-force)
- Retrieving the max amount of articles immediately, or rather a top #N articles. 
- Consider downloading all pmc articles
- Move database somewhere else?
- look into [graph tool](https://graph-tool.skewed.de/)
    - Hierarchical Stochastic Blockmodel Inference
    - Maximum Flow
- query for finding (abstract) attributes of articles related to sometihng
- query for finding relations between different labels connected by article
- figure out methods similairty score calc, barrelfield, barrel cortetx-- same thing (semantics)
  in our list we have barrel field, but people will search for barrel cortex
  similarity score is about keywords . When words SEEM similar, collapse them (but, then, still make it checkable..)
- add out/in/log/error dirs under /data/ and adjust scripts accordingly

### CONNECT GEPHI
    https://neo4j.com/labs/apoc/4.1/export/gephi/

### Bio.Medline.Record issue
Should we keep the mnemonics as attribute description in the database or translate them all to their corresponding interpretations, for example `MH` means `MeSH Terms`, which one do we use?


### (Neo4j) pitfalls
- don't use nodeId, dont store nodeId externally!

### notes for later (for creating package)
maybe watch [this](https://youtu.be/GIF3LaRqgXo?t=525)
- learn how to use pytest for testing
##### a setup.py file
probably  just use pycharm to make this
```python 
from setuptools import setup

setup(
    name='packagename',
    version='0.0.1',
    description='description of package',
    py_modules=['module1.py1', 'module2.py'],
    package_dir={'':'braintaxmap'}
)
```
- `python setup.py bdist_wheel`
- `pip install -e .`  in dir with package (`setup.py` file)