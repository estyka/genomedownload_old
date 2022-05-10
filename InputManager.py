import os
import os.path
import SharedConsts as sc
import pickle
from utils import logger
from Bio import Entrez


class InputManager:

    def __clean_input(self, bacteria_input):
        result = bacteria_input.lower().replace("_", " ")
        return result

    def validate_bacteria_input(self, email_address, bacteria_input):
        bacteria_input = self.__clean_input(bacteria_input)
        logger.info(f'validating organism input')
        return len(self.__get_assembly_ids(email_address, bacteria_input)) > 0

    def __get_assembly_ids(self, email_address, term):
        Entrez.email = email_address
        handle = Entrez.esearch(db="assembly", term=term, retmax='10000')
        record = Entrez.read(handle)
        ids = record['IdList']
        return ids

