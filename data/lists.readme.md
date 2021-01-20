# `DSM-5`

removed:
- All headings/subheadings
- Appendix

synonyms that stood between paranthesis are on a new line
- Something like "Schizotypal (Personality) Disorder" is broken up into 'Schizotypal  Disorder'  and 'Schizotypal Personality Disorder'
- entries with a slash are separated into two things Substance/Medication-Induced Mental Disorders -> 
`Substance Mental Disorders` and
`Medication-Induced Mental Disorders`
- summaries like Unspecified Sedative-, Hypnotic-, or Anxiolytic-Related Disorder
```python
def handle(line):

    if "/" in line:
        items = line.split("/")
        special = slashreg.findall(line)[0]
        a, b = special.split('/')
        a1 = line.replace(f"{a}/", '')
        b1 = line.replace(f"/{b}", '')
        return [line,a1,b1]
        
    if "(" in line:
        if line.endswith(")"):
            items = line.split("(")
            item = items[1][:-1]
            return [items[0], item]
        else:
            item = line[line.find("(")+1:line.find(")")]
            without = line.replace(f"({item}) ", "")
            return [line, without]

    else:
        return [line]
```