from ftplib import FTP
import time, os
import os.path
import SharedConsts as sc
import pickle
from datetime import datetime

SERVER_PATH = r"ftp.ncbi.nlm.nih.gov"
SERVER_BACTERIA_LOCATION = r"/genomes/refseq/bacteria/"
NUMBER_OF_TRIES = 500


class InputManager:

    def __clean_input(self, bacteria_input):
        result = bacteria_input.lower().replace("_", " ")
        return result

    def __get_bacteria_list(self):
        start = time.time()

        for i in range(NUMBER_OF_TRIES):
            try:
                print(f'starting num_try {i}')

                ftp = FTP(SERVER_PATH)
                ftp.login()
                print('logged in')

                ftp.cwd(SERVER_BACTERIA_LOCATION)
                print(f'changed dir to {SERVER_BACTERIA_LOCATION}')

                all_ncbi_bacterias = ftp.nlst(SERVER_BACTERIA_LOCATION)
                ftp.quit()  # This is the “polite” way to close a connection
                print("Successfully quited ftp connection")
                print(f"Runtime is {(time.time() - start)} seconds")
                return all_ncbi_bacterias

            except Exception as e:
                print(e)
                print(f'finished try number: {i} out of {NUMBER_OF_TRIES}')

    def __need_to_update_bacteria_list(self):
        bacteria_list_path = sc.PATH2BACTERIAS_LIST
        if not os.path.isfile(sc.PATH2BACTERIAS_LIST): #file doesn't exist yet
            return True
        last_modified_date = datetime.fromtimestamp(os.path.getmtime(bacteria_list_path)) #make sure this is created date if the file wasn't ever modified
        now = datetime.now()
        delta = now - last_modified_date
        if delta.days > 0:  # hasn't been updated today
            return True
        return False

    def validate_bacteria_input(self, bacteria_input):
        bacteria_list_path = sc.PATH2BACTERIAS_LIST
        if self.__need_to_update_bacteria_list():
            bacteria_list = self.__get_bacteria_list()
            with open(bacteria_list_path, 'wb') as f:
                pickle.dump(bacteria_list, f)
        else:  # no need to update list
            with open(bacteria_list_path) as f:
                bacteria_list = pickle.load(f)
        bacterias_of_interest = self.__get_bacteria_of_interest_paths_list(bacteria_input, bacteria_list)
        return len(bacterias_of_interest) > 0

    def __get_bacteria_of_interest_paths_list(self, bacteria_input, bacteria_list):
        bacterias_of_interest = [bacteria_path for bacteria_path in bacteria_list if
                                 self.__clean_input(os.path.basename(bacteria_path)).startswith(bacteria_input)]
        # TODO: make sure it's correct to do startswith
        return bacterias_of_interest

