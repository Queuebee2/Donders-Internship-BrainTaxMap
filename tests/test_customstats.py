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


stats = StatTracker()

stats.add(**{'val1':5, 'val2':6})
stats.add(bob=5, zob=54)
stats.add(bob=100)
test_dict = {'somethinsomething':554, "anotherthing":45}
stats.add(**test_dict)
print(stats)