import spacy
import textacy 
import nltk 

import time
ABSTRACT = """Toll-like receptor 4 (TLR4) is a crucial receptor in the innate immune system, and increasing evidence supports its role in inflammation, stress, and tissue injury, including injury to the lung and brain. We aimed to investigate the effects of TLR4 on neuroinflammation due to the lung-brain interaction in mechanically ventilated mice. Male wild-type (WT) C57BL/6 and TLR4 knockout (TLR4 KO) mice were divided into three groups: (1) control group (C): spontaneous breathing; (2) anesthesia group (A): spontaneous breathing under anesthesia; and (3) mechanical ventilation group (MV): 6h of MV under anesthesia. The behavioral responses of mice were tested with fear conditioning tests. The histological changes in the lung and brain were assessed using hematoxylin-eosin (HE) staining. The level of TLR4 mRNA 
in tissue was measured using reverse transcription-polymerase chain reaction (RT-PCR). The levels of inflammatory cytokines were measured with an enzyme-linked immunosorbent assay (ELISA). Microgliosis, astrocytosis, and the TLR4 immunoreactivity in the hippocampus were measured by double immunofluorescence. MV mice exhibited 
impaired cognition, and this impairment was less severe in TLR4 KO mice than in WT mice. In WT mice, MV increased TLR4 mRNA expression in the lung and brain. MV induced mild lung injury, which was prevented in TLR4 KO mice. MV mice exhibited increased levels of inflammatory cytokines, increased microglia and astrocyte activation. Microgliosis was alleviated in TLR4 KO mice. MV mice exhibited increased TLR4 immunoreactivity, which was expressed in microglia and astrocytes. These results demonstrate that TLR4 is involved in neuroinflammation due to the lung-brain interaction and that TLR4 KO loves neuroinflammation due to lung-brain interaction after prolonged MV. In addition, Administration of a TLR4 antagonist (100mug/mice) to WT mice also significantly attenuated neuroinflammation of lung-brain interaction due to prolonged MV. TLR4 antagonism may be a new and novel approach for the treatment and management of neuroinflammation in long-term mechanically ventilated patients.
"""

"""currently testing:
 the word 'ameliorates' was not found with the simple language model, maybe the large one (750mB!!!) 
 understands it
"""
tstart = time.perf_counter()


nlp = spacy.load("en_core_web_sm")
for sentence in nltk.sent_tokenize(ABSTRACT):
    doc = nlp(sentence)
    svo_triples = textacy.extract.subject_verb_object_triples(doc)
    svo_triples = list(svo_triples)
    if len(svo_triples) > 0:
        print(svo_triples)
tend = time.perf_counter()
duration = tend-tstart
print(f'Small model took {duration} seconds')


tstart = time.perf_counter()
#  python -m spacy download en_core_web_lg 750Mb
nlp = spacy.load("en_core_web_lg")
for sentence in nltk.sent_tokenize(ABSTRACT):
    doc = nlp(sentence)
    svo_triples = textacy.extract.subject_verb_object_triples(doc)
    svo_triples = list(svo_triples)
    if len(svo_triples) > 0:
        print(svo_triples)
tend = time.perf_counter()

duration = tend-tstart
print(f'Large model took {duration} seconds')
