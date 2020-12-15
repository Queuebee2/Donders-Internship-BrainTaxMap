
from dotenv import load_dotenv, find_dotenv
import dotenv

load_dotenv(find_dotenv())

dotenv.get_key(find_dotenv(), 'KEYTESTER')

dotenv.set_key(find_dotenv(), 'KEYTESTER', '.ENV_VALUE_TEST')

print('keytester value:' , dotenv.get_key(find_dotenv(), 'KEYTESTER'))

dotenv.unset_key(find_dotenv(), 'KEYTESTER')


dotenv.get_key(find_dotenv(), 'KEYTESTER')