import os
from pathlib import Path
from utils import State
from enum import Enum

# OUTPUT consts
K_MER_COUNTER_MATRIX_FILE_NAME = Path('CounterMatrixForUI.csv')
RESULTS_FOR_OUTPUT_CLASSIFIED_RAW_FILE_NAME = Path('ResultsForPostProcessClassifiedRaw.csv')
RESULTS_SUMMARY_FILE_NAME = Path('summary_results.txt')
RESULTS_FOR_OUTPUT_UNCLASSIFIED_RAW_FILE_NAME = Path('ResultsForPostProcessUnClassifiedRaw.csv')
TEMP_CLASSIFIED_IDS = Path('TempClassifiedIds.txt')
TEMP_UNCLASSIFIED_IDS = Path('TempUnClassifiedIds.txt')
INPUT_CLASSIFIED_FILE_NAME = Path('classified.fasta')
INPUT_UNCLASSIFIED_FILE_NAME = Path('unclassified.fasta')
FINAL_OUTPUT_DIR_NAME = Path(os.path.join("results", "fasta.zip"))
#FINAL_OUTPUT_FILE_NAME = Path('FilteredResults.txt.gz')
KRAKEN_SUMMARY_RESULTS_FOR_UI_FILE_NAME = Path('summary_stat_UI.json')
RANK_KRAKEN_TRANSLATIONS = {'U': 'Unclassified', 'R': 'Root', 'D': 'Domain', 'K': 'Kingdom', 'P': 'Phylum',
                            'C': 'Class', 'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'}

PATH_TO_OUTPUT_PROCESSOR_SCRIPT = Path("/groups/pupko/alburquerque/NgsReadClearEngine/OutputProcessor.py")
DF_LOADER_CHUCK_SIZE = 1e6
RESULTS_COLUMNS_TO_KEEP = ['is_classified', 'read_name', 'max_specie', 'classified_species', 'read_length', 'max_k_mer_p',
                           'all_classified_K_mers', 'split']
SUMMARY_RESULTS_COLUMN_NAMES = ['percentage_of_reads', 'number_of_reads_under', 'number_of_reads_exact', 'rank_code',
                                'ncbi_taxonomyID', 'name']
UNCLASSIFIED_COLUMN_NAME = 'Non Bacterial'
KRAKEN_UNCLASSIFIED_COLUMN_NAME = 'unclassified (taxid 0)'
UNCLASSIFIED_BACTERIA_NAME = 'Unclassified Bacteria'
UNCLASSIFIED_BACTERIA_ID = '-1'

# PBS Listener consts
JOB_NUMBER_COL = 'job_number'
JOB_NAME_COL = 'job_name'
JOB_STATUS_COL = 'job_status'
JOB_ELAPSED_TIME = 'elapsed_time'
JOB_CHANGE_COLS = [JOB_NUMBER_COL, JOB_NAME_COL, JOB_STATUS_COL]
QstatDataColumns = [JOB_NUMBER_COL, 'username', 'queue', JOB_NAME_COL, 'session_id', 'nodes', 'cpus', 'req_mem',
                    'req_time', JOB_STATUS_COL, JOB_ELAPSED_TIME]
SRVER_USERNAME = 'bioseq'
JOB_RUNNING_TIME_LIMIT_IN_HOURS = 10

# Job listener and management function naming
LONG_RUNNING_JOBS_NAME = 'LongRunning'
QUEUE_JOBS_NAME = 'Queue'
NEW_RUNNING_JOBS_NAME = 'NewRunning'
FINISHED_JOBS_NAME = 'Finished'
ERROR_JOBS_NAME = 'Error'
WEIRD_BEHAVIOR_JOB_TO_CHECK = ''
PATH2SAVE_PROCESS_DICT = r'/data/www/flask/genomedownload/SavedObjects/processes.dict'
PATH2BACTERIAS_LIST = r'/data/www/flask/genomedownload/SavedObjects/ncbi_bacterias.list'

INTERVAL_BETWEEN_LISTENER_SAMPLES = 5  # in seconds

PATH_2_DOWNLOAD_SCRIPT = r"/groups/pupko/estykatzeff/scripts/download_files_from_ncbi.py"
#PATH_2_DOWNLOAD_SCRIPT = r"/data/www/flask/genomedownload/PBS_process/download_files_from_ncbi.py"

JOB_PREFIX = 'GD'
POSTPROCESS_JOB_PREFIX = 'PP'
#KRAKEN_JOB_PREFIX = ''

# todo: replace the conda env
GENOME_DOWNLOAD_PROCESS_TEMPLATE = '''
#!/bin/bash -x
#PBS -S /bin/bash
#PBS -q lifesciweb
#PBS -N {job_name}
#PBS -e {error_files_path}
#PBS -o {output_files_path}
#PBS -r y
hostname
echo job_name: {job_name}
echo $PBS_JOBID
module list
module load python/python-anaconda3.6.5;python {path_to_python_script} {path_results} {speices_to_download}
'''

# post processing

class EMAIL_CONSTS:
    FINISHED_TITLE = f'Genome Downloader - Job finished'
    FINISHED_CONTENT = '''Thanks, for using Genome Downloader\nYour results are at:\nhttp://genomedownload.tau.ac.il/process_state/{process_id}\nPlease, remember to cite us'''
    CRASHED_TITLE = f'Genome Downloader - Job crashed'
    CRASHED_CONTENT = '''Thanks, for using Genome Downloader\nYour results are at:\nhttp://genomedownload.tau.ac.il/process_state/{process_id}\nPlease, remember to cite us'''




class UI_CONSTS:
    static_folder_path = 'gifs/'
    
    states_gifs_dict = {
    State.Running: {
        "background": "#9db09f",
        "gif_id": "aiqIqtW2utnkk"
    },
    State.Finished: {
        "background": "#1674d2",
        "gif_id": "TvLuZ00OIADoQ"
    },
    State.Crashed: {
        "background": "#1674d2",
        "gif_id": "TvLuZ00OIADoQ"
    },
    State.Waiting:  {
        "background": "#1674d2",
        "gif_id": "TvLuZ00OIADoQ"
    },
    State.Init:  {
        "background": "#1674d2",
        "gif_id": "TvLuZ00OIADoQ"
    },
    State.Queue:  {
        "background": "#1674d2",
        "gif_id": "TvLuZ00OIADoQ"
    },
    }
    
    states_text_dict = {
        State.Running: "Your process is running",
        State.Finished: "Your process finished... Redirecting to results page", #TODO is needed??
        State.Crashed: "Your process crashed\n we suggest you rerun the process.", #TODO finish
        State.Waiting: "We currently run other processes :( \n Your process will start soon",
        State.Init: "We are verifing your input, your process will start shortly",
        State.Queue: "Job is queued",
    }
    
    global allowed_files_str
    ALLOWED_EXTENSIONS = {'fasta', 'fastqc', 'gz'}
    allowed_files_str = ', '.join(ALLOWED_EXTENSIONS) #better to path string than list


    class UI_Errors(Enum):
        UNKNOWN_PROCESS_ID = 'The provided process id does not exist'
        INVALID_EXPORT_PARAMS ='invalid paramters for export'
        POSTPROCESS_CRASH = 'can\'t postprocess'
        INVALID_MAIL = 'invalid mail'
        CANT_ADD_PROCESS = 'can\'t add search process'
        INVALID_FILE = f'invalid file or file extenstion, please use a valid: {allowed_files_str} file'
        EXPORT_FILE_UNAVAILABLE = f'failed to export file, try to rerun the file'
        PAGE_NOT_FOUND = 'The requested page does not exist'

    PROCESS_INFO_PP = "We are processing your request, This may take several minutes. You may close this window, An email will be sent upon completion"
    PROCESS_INFO_KR = "We are processing your request, This may take several minutes for small files and several hours for larger ones. Please close this window, An email will be sent upon completion"