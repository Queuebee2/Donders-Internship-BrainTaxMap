
for item in [1,2,3,4,5]:

    lots_of_errors_loop = iter(range(item*10))
    
    c=0

    while True:
        try:
            item = next(lots_of_errors_loop)
            print(item)
            if item % 3 == 0:
                raise FileNotFoundError
            c+=1
            # do inner loop shit
        except StopIteration:
            print('stopiteration')
            break
        except Exception as e:
            print('setohming')
    

    