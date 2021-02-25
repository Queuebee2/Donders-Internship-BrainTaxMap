# Current branch: Internship end state
This is the final state of the repository at the end of my internship. 

# BrainTaxMap
The project goal is to map the behavioural taxonomy of the brain by systematically searching and analysing biomedical texts. At this stage, we explored retrieving all articles related to all neural structure names from the [Mouse Brain Atlas Ontology's StructureGraph](http://api.brain-map.org/api/v2/structure_graph_download/10.json) and then creating histograms with counts of the relations found with regular expressions, looking for combinations of structure names (in [`data/lists/structures`](https://github.com/Queuebee2/Donders-Internship-BrainTaxMap/blob/Internship-endstate/data/lists/structures)), verbs (in [`../verbs`](https://github.com/Queuebee2/Donders-Internship-BrainTaxMap/blob/Internship-endstate/data/lists/verbs)) and disorders (in [`../disorders`](https://github.com/Queuebee2/Donders-Internship-BrainTaxMap/blob/Internship-endstate/data/lists/disorders)) in sentences for each available abstract. Disorders from [ICD-11](https://icd.who.int/browse11/l-m/en) are used.



