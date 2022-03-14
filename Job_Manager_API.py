import os
import shutil
import uuid
import json
import zipfile
from Job_Manager_Thread_Genome_Download import Job_Manager_Thread_Genome_Download
from utils import send_email, State, logger, LOGGER_LEVEL_JOB_MANAGE_API
from SharedConsts import EMAIL_CONSTS, FINAL_OUTPUT_DIR_NAME
from InputManager import InputManager
logger.setLevel(LOGGER_LEVEL_JOB_MANAGE_API)

class Job_Manager_API:
    def __init__(self, max_number_of_process: int, upload_root_path: str, func2update_html):
        self.__upload_root_path = upload_root_path
        self.__j_manager = Job_Manager_Thread_Genome_Download(max_number_of_process, upload_root_path, self.__process_state_changed)
        self.__func2update_html = func2update_html
        self.input_manager = InputManager()

    def __build_and_send_mail(self, process_id, subject, content, email_address):
        try:
            send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>',
                       email_address, subject=subject,
                       content= content)
            logger.info(f'sent email to {email_address}')
        except:
            logger.exception(f'failed to sent email to {email_address}')

    def __process_state_changed(self, process_id, state, email_address):
        if state == State.Finished:
            if email_address != None:
                self.__build_and_send_mail(process_id, EMAIL_CONSTS.FINISHED_TITLE, EMAIL_CONSTS.FINISHED_CONTENT.format(process_id=process_id), email_address)
        elif state == State.Crashed:
            if email_address != None:
                self.__build_and_send_mail(process_id, EMAIL_CONSTS.CRASHED_TITLE, EMAIL_CONSTS.CRASHED_CONTENT.format(process_id=process_id), email_address)
        self.__func2update_html(process_id, state)

    def __delete_folder(self, process_id):
        logger.info(f'process_id = {process_id}')
        folder2remove = os.path.join(self.__upload_root_path, process_id)
        shutil.rmtree(folder2remove)

    def __validate_organism_name(self, process_id, organism_name):
        parent_folder = os.path.join(self.__upload_root_path, process_id)
        if not os.path.isdir(parent_folder):
            logger.warning(f'process_id = {process_id} doen\'t have a dir')
            return False
        if self.input_manager.validate_bacteria_input(organism_name):
            return True
        self.__delete_folder(process_id)
        logger.warning(f'validation failed {organism_name}, deleting process folder')
        return False
        
    def __validate_email_address(self, email_address):
        if len(email_address) > 100:
            return False
        if '@' in email_address and '.' in email_address:
            return True
        return False

    def get_new_process_id(self):
        return str(uuid.uuid4())

    def add_process(self, process_id: str, email_address: str, organism_name: str):
        logger.info(f'process_id = {process_id} email_address = {email_address} organism_name = {organism_name}')
        is_valid_organism = self.__validate_organism_name(process_id, organism_name)
        is_valid_email = self.__validate_email_address(email_address)
        if is_valid_organism and is_valid_email:
            logger.info(f'validated organism name and email address')
            self.__j_manager.add_download_process(process_id, email_address, organism_name)
            return True
        logger.warning(f'process_id = {process_id}, can\'t add process: is_valid_email = {is_valid_email}')
        return False
            
    def export_dir(self, process_id: str):
        parent_folder = os.path.join(self.__upload_root_path, process_id)
        if os.path.isdir(parent_folder):
            dir2return = os.path.join(parent_folder, FINAL_OUTPUT_DIR_NAME)
            z = zipfile.ZipFile(dir2return)
            logger.info(f'z.namelist: {z.namelist()}')
            if z.namelist() is not None:
                return dir2return
        logger.warning(f'process_id = {process_id} doen\'t have a result dir2return: {dir2return}')
        return None
    
    def get_kraken_job_state(self, process_id):
        return self.__j_manager.get_job_state(process_id)
        
        
    def get_UI_matrix(self, process_id):
        parent_folder = os.path.join(self.__upload_root_path, process_id)
        csv_UI_matrix_path = os.path.join(parent_folder, K_MER_COUNTER_MATRIX_FILE_NAME)
        summary_stats_json_path = os.path.join(parent_folder, KRAKEN_SUMMARY_RESULTS_FOR_UI_FILE_NAME)
        df2return = None
        json2return = None
        if os.path.isfile(csv_UI_matrix_path):
            df2return = pd.read_csv(csv_UI_matrix_path,index_col=0)
        if os.path.isfile(summary_stats_json_path):
            json2return = json.load(open(summary_stats_json_path))
        
        logger.info(f'process_id = {process_id} df2return = {df2return} json2return = {json2return}')
        return df2return, json2return
