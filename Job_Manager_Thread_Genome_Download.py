import os
import SharedConsts as sc
from PBS_process.dowload_process import create_download_process
from Job_Manager_Thread_Safe import Job_Manager_Thread_Safe
from utils import logger


class Job_Manager_Thread_Genome_Download:
    def __init__(self, max_number_of_process: int, upload_root_path: str, func2update_html_download):
        self.__func2update_html_download = func2update_html_download
        function2call_processes_changes_state = {
            sc.JOB_PREFIX: self.__func2update_html_download,
        }
        function2append_process = {
            sc.JOB_PREFIX: self.__download_process,
        }
        paths2verify_process_ends = {
            #when the job crashes/ finished this file path will be checked to set the change to finished if file exists of crashed if file doesn't.
            #for a string of: '' it won't set the state
            sc.JOB_PREFIX: lambda process_id: os.path.join(os.path.join(upload_root_path, process_id), sc.K_MER_COUNTER_MATRIX_FILE_NAME),
        }
        self.__job_manager = Job_Manager_Thread_Safe(max_number_of_process, upload_root_path, function2call_processes_changes_state, function2append_process, paths2verify_process_ends)

    def __download_process(self, process_folder_path: str, email_address: str, organism_name: str):
        logger.info(f'process_folder_path = {process_folder_path} email_address = {email_address}  organism_name = {organism_name}')
        pbs_id = create_download_process(process_folder_path, organism_name)
        print(f'pbs_id = {pbs_id}')
        return pbs_id
            
    def __get_state(self, process_id, job_prefix):
        state = self.__job_manager.get_job_state(process_id, job_prefix)
        if state:
            return state
        logger.warning(f'process_id = {process_id}, job_prefix = {job_prefix} not in __job_manager')
        return None
        
    def get_download_job_state(self, process_id):
        return self.__get_state(process_id, sc.JOB_PREFIX)
        
    
    def add_download_process(self, process_id: str, email_address: str, organism_name: str):
        logger.info(f'process_id = {process_id}')
        self.__job_manager.add_process(process_id, sc.JOB_PREFIX, email_address, organism_name)
    
    
    def get_job_state(self, process_id: str, job_prefix: str):
        logger.info(f'process_id = {process_id} job_prefix = {job_prefix}')
        return self.__job_manager.get_job_state(process_id, job_prefix)
