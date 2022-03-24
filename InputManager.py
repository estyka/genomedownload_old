import os
import os.path
import SharedConsts as sc
import pickle
from utils import logger


class InputManager:

    def __clean_input(self, bacteria_input):
        result = bacteria_input.lower().replace("_", " ")
        return result

    def validate_bacteria_input(self, bacteria_input):
        bacteria_input = self.__clean_input(bacteria_input)
        logger.info(f'validating organism input')
        bacteria_list_path = sc.PATH2BACTERIAS_LIST
        with open(bacteria_list_path, 'rb') as f:
            bacteria_list = pickle.load(f)
        bacterias_of_interest = self.__get_bacteria_of_interest_paths_list(bacteria_input, bacteria_list)
        logger.info(f'bacterias_of_interest = {bacterias_of_interest}')
        return len(bacterias_of_interest) > 0

    def __get_bacteria_of_interest_paths_list(self, bacteria_input, bacteria_list):
        bacterias_of_interest = [bacteria_path for bacteria_path in bacteria_list if
                                 self.__clean_input(os.path.basename(bacteria_path)).startswith(bacteria_input)]
        # TODO: make sure it's correct to do startswith
        return bacterias_of_interest