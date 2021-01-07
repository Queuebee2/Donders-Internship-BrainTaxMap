from collections import defaultdict, Counter

bob = defaultdict(Counter)
for thing in 'qwertyuiop':
    for otherthing in 'azertypoi':
        bob[thing][otherthing] +=1

print(bob)

