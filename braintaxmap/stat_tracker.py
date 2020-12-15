class StatTracker():
    def __init__(self):
        self.statnames = []
        self.statvalues = []
        
    def add(self, **kwargs):
        for k, v in kwargs.items():
            self.statnames.append(k)
            self.statvalues.append(v)

    def __str__(self):
        for stat, val in zip(self.statnames, self.statvalues):
            print(f'{stat} | {val}')
        return '- - - - - - - - - - -'
