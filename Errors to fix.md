-----------------------------------------------------------------
http.client.IncompleteRead: IncompleteRead(0 bytes read)
trying 
https://stackoverflow.com/questions/41529016/python-http-client-incomplete-read0-bytes-read-error

adding this to 
`Bio.Medline.__init__.parse()`
 everything in a Try block
```python

    for line in handle:
        try:
            #....
        except IncompleteRead:
            print('encountered incomplete read error...')
            continue
```
-----------------------------------------------------------------
# error caused by pubmed, giving no id.
k id: v ['372060 Error occurred: The following PMID is not available: 33372060']
https://pubmed.ncbi.nlm.nih.gov/33372060/



-----------------------------------------------------------------
  File "c:\Users\Milai\Documents\Github\Donders-Internship-BrainTaxMap\braintaxmap\tools.py", line 49, in log_all_svos
    for svo, count in stats['allSVOs'].most_common():
KeyError: 'allSVOs'



Traceback (most recent call last):
  File "C:\Users\Milai\AppData\Local\Programs\Python\Python37\lib\logging\__init__.py", line 1025, in emit
    msg = self.format(record)
  File "C:\Users\Milai\AppData\Local\Programs\Python\Python37\lib\logging\__init__.py", line 869, in format
    return fmt.format(record)
  File "C:\Users\Milai\AppData\Local\Programs\Python\Python37\lib\logging\__init__.py", line 608, in format
    record.message = record.getMessage()
  File "C:\Users\Milai\AppData\Local\Programs\Python\Python37\lib\logging\__init__.py", line 369, in getMessage
    msg = msg % self.args
TypeError: not all arguments converted during string formatting
Call stack:
  File "c:\Users\Milai\Documents\Github\Donders-Internship-BrainTaxMap\main.py", line 433, in <module>
    logger.exception(f'error for query "{query}"" at rec_index "{rec_index}"', e)
Message: 'error for query "barrel cortex"" at rec_index "1507"'
Arguments: (ConnectionResetError(10054, 'An existing connection was forcibly closed by the remote host', None, 10054, None),)

-----------------------------------------------------------------

RuntimeError: Search Backend failed: read request has timed out. peer: 