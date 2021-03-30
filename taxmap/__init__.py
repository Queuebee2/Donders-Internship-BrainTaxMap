import logging
import os
from pathlib import Path

from dotenv import find_dotenv, load_dotenv

logger = logging.getLogger(__name__)

# load environment variables from .env file
load_dotenv(find_dotenv())
DEV_EMAIL = os.environ.get('DEV_EMAIL') or None
NCBI_API_KEY = os.getenv('NCBI_API_KEY') or None

# create paths
Path('data/').mkdir(parents=True, exist_ok=True)
Path('data/output/graphs/treemaps/').mkdir(parents=True, exist_ok=True)
Path('data/output/records/').mkdir(parents=True, exist_ok=True)

DATA_DIR = Path('data/').resolve()
DEFAULT_RECORDS_DIR = Path('data/output/records/').resolve()

from .const import *
