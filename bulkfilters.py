
import abc

class Filter(metaclass=abc.ABCMeta):
    def __name__(self):
        return self.__class__.__name__
    
    @abc.abstractmethod
    def filter(self, *args, **kwargs):
        pass

class RecordFilter(Filter):
    def __init__(self, record, state=1, *args, **kwargs):
        self.value = self.filter(record, *args, **kwargs) if state else 0

class NoAbstract(RecordFilter):
    def filter(self, record):
        return 'AB' not in  record.keys()

class HasRelevantMeshterms(RecordFilter):

    def filter(self, record, includelist, excludelist): # or just one record?
        """ find excluded/included mesh terms in the meshterm list
        and return those lists
        """
        included, excluded = list(), list()
        if 'MH' in record.keys():
            for m in record['MH']:
                meshterm = m.lower()
                for incl in includelist:
                    if incl.lower() in meshterm:
                        included.append(meshterm)
                
                for excl in excludelist:
                    if excl.lower() in meshterm:
                        excluded.append(meshterm)

        return included, excluded

class RelevantAbstract(RecordFilter):

    def filter(self, record, sets_with_relevant_terms, force_one_of_each=True):
        """within a string of text, checks whether a word from a list/set exists in it.
        can be given multiple lists. If the force_one_of_each parameter is True,
        one word from each list has to be included in the text to return True
        """
        if 'AB' not in record.keys():
            return None
        
        checks = [0 for _ in sets_with_relevant_terms]

        for word in record['AB'].split(' '):
            for i, s in enumerate(sets_with_relevant_terms):
                for keyword in s:
                    if force_one_of_each:
                        checks[i]=1
                        break # this set is now included, skip to next
                    else:
                        return True
        else:
            if force_one_of_each and all(checks):
                return True
            else:
                return False
