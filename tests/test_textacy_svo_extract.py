import spacy
import textacy 
import nltk 
import time

import time
ABSTRACT = """Toll-like receptor 4 (TLR4) is a crucial receptor in the innate immune system, and increasing evidence supports its role in inflammation, stress, and tissue injury, including injury to the lung and brain. We aimed to investigate the effects of TLR4 on neuroinflammation due to the lung-brain interaction in mechanically ventilated mice. Male wild-type (WT) C57BL/6 and TLR4 knockout (TLR4 KO) mice were divided into three groups: (1) control group (C): spontaneous breathing; (2) anesthesia group (A): spontaneous breathing under anesthesia; and (3) mechanical ventilation group (MV): 6h of MV under anesthesia. The behavioral responses of mice were tested with fear conditioning tests. The histological changes in the lung and brain were assessed using hematoxylin-eosin (HE) staining. The level of TLR4 mRNA 
in tissue was measured using reverse transcription-polymerase chain reaction (RT-PCR). The levels of inflammatory cytokines were measured with an enzyme-linked immunosorbent assay (ELISA). Microgliosis, astrocytosis, and the TLR4 immunoreactivity in the hippocampus were measured by double immunofluorescence. MV mice exhibited 
impaired cognition, and this impairment was less severe in TLR4 KO mice than in WT mice. In WT mice, MV increased TLR4 mRNA expression in the lung and brain. MV induced mild lung injury, which was prevented in TLR4 KO mice. MV mice exhibited increased levels of inflammatory cytokines, increased microglia and astrocyte activation. Microgliosis was alleviated in TLR4 KO mice. MV mice exhibited increased TLR4 immunoreactivity, which was expressed in microglia and astrocytes. These results demonstrate that TLR4 is involved in neuroinflammation due to the lung-brain interaction and that TLR4 KO ameliorates neuroinflammation due to lung-brain interaction after prolonged MV. In addition, Administration of a TLR4 antagonist (100mug/mice) to WT mice also significantly attenuated neuroinflammation of lung-brain interaction due to prolonged MV. TLR4 antagonism may be a new and novel approach for the treatment and management of neuroinflammation in long-term mechanically ventilated patients.
"""

"""currently testing:
 the word 'ameliorates' was not found with the simple language model, maybe the large one (750mB!!!) 
 understands it
"""
tstart = time.perf_counter()

AB2= """PURPOSE/AIM OF THE STUDY: Whiskers are important sensory organs that play a key role in rodents' discriminative and exploration behaviours and unilateral injuries of the somatosensory cortex related to whisker barrel cortex can change the activity of neurons in the intact contralateral barrel cortex. We evaluated the effects of unilateral mechanical lesion of right barrel cortex on novel texture discrimination in behavioural test and neuronal responses of left barrel cortex. MATERIALS AND METHODS: Ten days after a unilateral mechanical lesion in the right barrel cortex, adult male rats were experimented regarding three paired different textures in novel texture discrimination test dependent on whiskers. In addition, responses of left barrel cortical neurons to controlled deflections of right whiskers were recorded using extracellular single-unit recordings technique. RESULTS: Data analysis showed that the discrimination ratio and preference indexes as criteria to find a novel texture significantly decreased in the lesion group compared to the intact rats (p < .05). In electrophysiological level, the barrel neural cortical spontaneous activity and the ON and OFF response magnitude of intact barrel cortex neurons in the lesion group decreased compared to the intact group (p < .05). CONCLUSIONS: The present study showed that unilateral mechanical lesion in the rats' barrel cortex cause a decrease in their abilities for discriminating textures, as well as, the anaesthetized rats whose response properties of intact barrel cortical area changed to whisker deflection, too. These changes can influence on the ability of rats to differentiate textures.
"""
nlp_large = spacy.load("en_core_web_lg") #  python -m spacy download en_core_web_lg 750Mb
nlp_small = spacy.load("en_core_web_sm")

for sent in nltk.sent_tokenize(AB2):
    print(sent)
    t = nlp_large(sent)
    t = list(textacy.extract.subject_verb_object_triples(t))
    print('larg',t)
    t = nlp_small(sent)
    t = list(textacy.extract.subject_verb_object_triples(t))
    print('sml',t)

def test_parsing_timedifference():
    AMOUNT=1000
    total_small=0
    total_large=0

    for i in range(1000):
        tstart = time.perf_counter()
        for sentence in nltk.sent_tokenize(ABSTRACT):
            doc_small = nlp_small(sentence)

            svo_triples_small = list(textacy.extract.subject_verb_object_triples(doc_small))

        tend = time.perf_counter()
        duration = tend-tstart
        # print(f'small model duration: {duration}')
        total_small+=duration  
            



        tstart = time.perf_counter()

        for sentence in nltk.sent_tokenize(ABSTRACT):

            doc_large = nlp_large(sentence)

            svo_triples_large = list(textacy.extract.subject_verb_object_triples(doc_large))

        tend = time.perf_counter()
        duration = tend-tstart
        # print(f'large model duration: {duration}')
        total_large+=duration

    print('total small:', total_small, 'average', total_small/AMOUNT, 'per run')
    print('total large:', total_large, 'average', total_large/AMOUNT, 'per run')


def test_loadingtimes():
    AMOUNT=10
    total_small=0
    total_large=0

    for i in range(10):
        tstart = time.perf_counter()
        nlp_small = spacy.load("en_core_web_sm")
        tend = time.perf_counter()
        duration = tend-tstart
        total_small+=duration  

        tstart = time.perf_counter()
        nlp_large = spacy.load("en_core_web_lg") 
        tend = time.perf_counter()
        duration = tend-tstart
        total_large+=duration

    print('loading total small:', total_small/AMOUNT)
    print('loading total large:', total_large/AMOUNT)

# test_parsing_timedifference()
# test_loadingtimes()