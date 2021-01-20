import requests
import dotenv

# taken and adapted from the example here
# https://github.com/ICD-API/Python-samples/blob/master/sample.py

def fetch_token(client_id, client_secret):
    token_endpoint = 'https://icdaccessmanagement.who.int/connect/token'
    scope = 'icdapi_access'
    grant_type = 'client_credentials'

    # set data to post
    payload = {
        'client_id'    : client_id, 
        'client_secret': client_secret, 
        'scope'        : scope, 
        'grant_type'   : grant_type
        }
            
    # make request
    r = requests.post(token_endpoint, data=payload, verify=False).json()
    token = r['access_token']
    return token



if __name__ == '__main__':
    client_id = dotenv.get_key(dotenv.find_dotenv(), 'ICD_ID')
    client_secret = dotenv.get_key(dotenv.find_dotenv(), 'ICD_CLIENT_SEC')
    token = fetch_token(client_id, client_secret)

    # access ICD API
    uri = 'http://id.who.int/icd/entity/334423054'

    # HTTP header fields to set
    headers = {
        'Authorization'  : 'Bearer ' + token, 
        'Accept'         : 'application/json', 
        'Accept-Language': 'en',
        'API-Version'    : 'v2'
        }
        
    # make request     

    r = requests.get(uri, headers=headers, verify=False)

    # print the result
    print(r.text)
