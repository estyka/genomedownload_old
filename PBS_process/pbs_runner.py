import pathlib
import subprocess
from subprocess import PIPE
import os
from SharedConsts import GENOME_DOWNLOAD_PROCESS_TEMPLATE, JOB_PREFIX, PATH_2_DOWNLOAD_SCRIPT, PATH_2_DOWNLOAD_LIST_SCRIPT, SPECIES_LIST_DOWNLOAD_PROCESS_TEMPLATE, PATH2BACTERIAS_LIST
from utils import logger

def create_download_process(input_path, species2download):
    # create the job
    species2download = species2download.replace(" ", "_") # turning parameter into one word
    job_unique_id = str(pathlib.Path(input_path).stem)
    temp_script_path = pathlib.Path().resolve() / f'temp_downloadgenome_file_{job_unique_id}.sh'
    job_name = f'{JOB_PREFIX}_{job_unique_id}'
    job_logs_path = str(pathlib.Path(input_path)) + '/'
    temp_script_text = GENOME_DOWNLOAD_PROCESS_TEMPLATE.format(job_name=job_name, error_files_path=job_logs_path,
                                                           output_files_path=job_logs_path, path_to_python_script=PATH_2_DOWNLOAD_SCRIPT,
                                                           path_results=job_logs_path, speices_to_download=species2download)

    # run the job
    with open(temp_script_path, 'w+') as fp:
        fp.write(temp_script_text)
    logger.info(f'submitting job, temp_script_path = {temp_script_path}:')
    logger.debug(f'{temp_script_text}')
    terminal_cmd = f'/opt/pbs/bin/qsub {str(temp_script_path)}'
    job_run_output = subprocess.run(terminal_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    print(f'job_run_output = {job_run_output}')
    os.remove(temp_script_path)

    return job_run_output.stdout.decode('utf-8').split('.')[0]
    
def create_download_species_list_process(input_path):
    # create the job
    job_unique_id = str(pathlib.Path(input_path).stem)
    temp_script_path = pathlib.Path().resolve() / f'temp_download_species_list_file_{job_unique_id}.sh'
    job_name = f'{JOB_PREFIX}_{job_unique_id}'
    job_logs_path = str(pathlib.Path(input_path)) + '/'
    temp_script_text = SPECIES_LIST_DOWNLOAD_PROCESS_TEMPLATE.format(job_name=job_name, error_files_path=job_logs_path,
                                                           output_files_path=job_logs_path, path_to_python_script=PATH_2_DOWNLOAD_LIST_SCRIPT,
                                                           path_results=PATH2BACTERIAS_LIST)

    # run the job
    with open(temp_script_path, 'w+') as fp:
        fp.write(temp_script_text)
    logger.info(f'submitting job, temp_script_path = {temp_script_path}:')
    logger.debug(f'{temp_script_text}')
    terminal_cmd = f'/opt/pbs/bin/qsub {str(temp_script_path)}'
    job_run_output = subprocess.run(terminal_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    print(f'job_run_output = {job_run_output}')
    #os.remove(temp_script_path)
    
    return job_run_output.stdout.decode('utf-8').split('.')[0]