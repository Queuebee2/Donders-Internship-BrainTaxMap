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


