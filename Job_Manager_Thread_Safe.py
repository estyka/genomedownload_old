import datetime
import os
import pickle
from threading import Lock
from apscheduler.schedulers.background import BackgroundScheduler
from JobListener import PbsListener
import SharedConsts as sc
from utils import State, logger, LOGGER_LEVEL_JOB_MANAGE_THREAD_SAFE
import traceback
from PBS_process.pbs_runner import run_create_download_species_list_process

logger.setLevel(LOGGER_LEVEL_JOB_MANAGE_THREAD_SAFE)


class Job_State:
    def __init__(self, folder_path: str, jobs_prefixes_lst: list, email_address):
        self.__folder_path = folder_path
        self.__job_states_dict = {prefix: None for prefix in jobs_prefixes_lst}
        self.__time_added = datetime.datetime.now()
        self.__pbs_id_dict = {prefix: None for prefix in jobs_prefixes_lst}
        self.__email_address = email_address

    def set_job_state(self, new_state: State, job_prefix: str):
        if job_prefix in self.__job_states_dict:
            self.__job_states_dict[job_prefix] = new_state
        else:
            logger.error(f'job_prefix = {job_prefix} not in self.__job_states_dict = {self.__job_states_dict}')
        
    def get_job_state(self, job_prefix: str):
        if job_prefix in self.__job_states_dict:
            return self.__job_states_dict[job_prefix]
        else:
            logger.error(f'job_prefix = {job_prefix} not in self.__job_states_dict = {self.__job_states_dict}')
        
    def set_pbs_id(self, pbs_id: str, job_prefix: str):
        if job_prefix in self.__job_states_dict:
            self.__pbs_id_dict[job_prefix] = pbs_id
        else:
            logger.error(f'job_prefix = {job_prefix} not in self.__pbs_id_dict = {self.__pbs_id_dict}')
        
    def get_pbs_id(self, job_prefix: str):
        if job_prefix in self.__pbs_id_dict:
            return self.__pbs_id_dict[job_prefix]
        else:
            logger.error(f'job_prefix = {job_prefix} not in self.__pbs_id_dict = {self.__pbs_id_dict}')
        
    def get_email_address(self):
        return self.__email_address

class Job_Manager_Thread_Safe:
    def __init__(self, max_number_of_process: int, upload_root_path: str, function2call_processes_changes_state: dict, function2append_process: dict, paths2verify_process_ends: dict):
        self.__max_number_of_process = max_number_of_process
        self.__upload_root_path = upload_root_path
        self.__processes_state_dict = self.__read_processes_state_dict2file()
        self.__clean_processes_state_dict()
        self.__mutex_processes_state_dict = Lock()
        self.__mutex_processes_waiting_queue = Lock()
        self.__waiting_list = []
        assert len(function2call_processes_changes_state) == len(function2append_process) == len(paths2verify_process_ends), f'verify function2call_processes_changes_state, function2append_process and paths2verify_process_ends have all the required job_prefixes. Their len should be the same'
        assert function2call_processes_changes_state.keys() == function2append_process.keys() == paths2verify_process_ends.keys(), f'verify function2call_processes_changes_state, function2append_process and paths2verify_process_ends have the same keys. It should contain the job_prefixes'
        self.jobs_prefixes_lst = list(function2call_processes_changes_state.keys())
        self.__function2append_process = function2append_process
        self.__paths2verify_process_ends = paths2verify_process_ends
        function_to_call_listener = {}
        for job_prefix in function2call_processes_changes_state.keys():
            function_to_call_listener[job_prefix] = self.__make_function_dict4listener(lambda process_id, state, _job_prefix=job_prefix: self.__set_process_state(process_id, state, _job_prefix, function2call_processes_changes_state[_job_prefix]))
        # create listener on queue
        self.__listener = PbsListener(function_to_call_listener)
        self.__scheduler = BackgroundScheduler()
        self.__scheduler.add_job(self.__listener.run, 'interval', seconds=sc.INTERVAL_BETWEEN_LISTENER_SAMPLES)
        self.__scheduler.add_job(run_create_download_species_list_process, 'interval', seconds=sc.INTERVAL_BETWEEN_BACTERIA_LIST_UPDATERS)
        #self.__scheduler.add_job(run_create_download_species_list_process, 'interval', seconds=60*60)  #testing 1 hour
        logger.info(f'scheduler added run_create_download_species_list_process job')
        self.__scheduler.start()
        
    def __save_processes_state_dict2file(self):
        file_to_store = open(sc.PATH2SAVE_PROCESS_DICT, "wb")
        pickle.dump(self.__processes_state_dict, file_to_store)
        file_to_store.close()
        
    def __read_processes_state_dict2file(self):
        dict2return = {}
        if os.path.isfile(sc.PATH2SAVE_PROCESS_DICT):
            file_to_read = open(sc.PATH2SAVE_PROCESS_DICT, "rb")
            dict2return = pickle.load(file_to_read)
            file_to_read.close()
        logger.info(f'dict2return = {dict2return}')
        return dict2return

    def __clean_processes_state_dict(self):
        for process_id in list(self.__processes_state_dict):
            folder_path = os.path.join(self.__upload_root_path, process_id)
            if not os.path.isdir(folder_path):
                del self.__processes_state_dict[process_id]

    def __calc_num_running_processes(self):
        running_processes = 0
        self.__mutex_processes_state_dict.acquire()
        for process_id in self.__processes_state_dict:
            for job_prefix in self.jobs_prefixes_lst:
                if self.__processes_state_dict[process_id].get_job_state(job_prefix) in [State.Running, State.Queue, State.Init]:
                    running_processes += 1
        self.__mutex_processes_state_dict.release()
        return running_processes

    def __calc_process_id(self, pbs_id):
        clean_pbs_id = pbs_id.split('.')[0]
        process_id2return = None
        self.__mutex_processes_state_dict.acquire()
        for process_id in self.__processes_state_dict:
            for job_prefix in self.jobs_prefixes_lst:
                if clean_pbs_id == self.__processes_state_dict[process_id].get_pbs_id(job_prefix):
                    process_id2return = process_id
                    break
        self.__mutex_processes_state_dict.release()
        if not process_id2return:
            logger.warning(f'clean_pbs_id = {clean_pbs_id} not in __processes_state_dict')
        return process_id2return

    def __log_and_set_change(self, pbs_id, set_process_state_func, state):
        process_id = self.__calc_process_id(pbs_id)
        logger.info(f'pbs_id = {pbs_id}  process_id = {process_id} state is {state}')
        set_process_state_func(process_id, state)

    def __make_function_dict4listener(self, set_process_state):
        return {
            sc.LONG_RUNNING_JOBS_NAME: lambda x: self.__log_and_set_change(x, set_process_state, State.Running), #TODO handle -currently same behevior as running
            sc.NEW_RUNNING_JOBS_NAME: lambda x: self.__log_and_set_change(x, set_process_state, State.Running),
            sc.QUEUE_JOBS_NAME: lambda x: self.__log_and_set_change(x, set_process_state, State.Queue),
            sc.FINISHED_JOBS_NAME: lambda x: self.__log_and_set_change(x, set_process_state, State.Finished),
            sc.ERROR_JOBS_NAME: lambda x: self.__log_and_set_change(x, set_process_state, State.Crashed),
        }

    def __set_process_state(self, process_id, state, job_prefix, func2update):
        logger.info(f'process_id = {process_id}, job_prefix = {job_prefix} state = {state}')
        email_address = None
        
        if (state == State.Finished or state == State.Crashed) and process_id in self.__processes_state_dict:
            dir2check = self.__paths2verify_process_ends[job_prefix](process_id)
            if dir2check != '':  # if dir2check is '' don't change the state
                if os.path.isdir(dir2check):
                    state = State.Finished
                else:
                    state = State.Crashed
        self.__mutex_processes_state_dict.acquire()
        if process_id in self.__processes_state_dict:
            self.__processes_state_dict[process_id].set_job_state(state, job_prefix)
            email_address = self.__processes_state_dict[process_id].get_email_address()
        else:
            # TODO handle
            logger.warning(f'process_id {process_id} not in __processes_state_dict: {self.__processes_state_dict}')
        self.__mutex_processes_state_dict.release()
        
        # don't put inside the mutex area - the funciton acquire the mutex too
        if state == State.Finished or state == State.Crashed:
            self.__add_process_from_waiting_list()
        
        #update file with current process_state_dict
        self.__save_processes_state_dict2file()
        func2update(process_id, state, email_address)

    def __add_process_from_waiting_list(self):
        process2add, job_type, running_arguments = self.__pop_from_waiting_queue()
        if process2add:
            logger.debug(f'adding new process after processed finished process2add = {process2add} job_type = {job_type}')
            self.add_process(process2add, job_type, *running_arguments)

    def add_process(self, process_id: str, job_prefix, *args):
        logger.info(f'process_id = {process_id}, job_prefix = {job_prefix}, args = {args}')
        # don't put inside the mutex area - the funciton acquire the mutex too
        running_processes = self.__calc_num_running_processes()
        self.__mutex_processes_state_dict.acquire()
        if running_processes < self.__max_number_of_process:
            process_folder_path = os.path.join(self.__upload_root_path, process_id)
            if process_id not in self.__processes_state_dict:
                email_address = args[0]
                self.__processes_state_dict[process_id] = Job_State(process_folder_path, self.__function2append_process.keys(), email_address)
            
            try:
                pbs_id = self.__function2append_process[job_prefix](process_folder_path, *args)
                logger.debug(f'process_id = {process_id} job_prefix = {job_prefix} pbs_id = {pbs_id}, process has started')
                self.__processes_state_dict[process_id].set_job_state(State.Init, job_prefix)
                self.__processes_state_dict[process_id].set_pbs_id(pbs_id, job_prefix)
            except Exception as e:
                traceback.print_exc()
                logger.error(e)

        else:
            logger.info(f'process_id = {process_id} job_prefix = {job_prefix}, adding to waiting list')
            self.__mutex_processes_waiting_queue.acquire()
            self.__waiting_list.append((process_id, job_prefix, args))
            self.__mutex_processes_waiting_queue.release()
            
        self.__mutex_processes_state_dict.release()
        #update file with current process_state_dict
        self.__save_processes_state_dict2file()
    
    def __pop_from_waiting_queue(self):
        self.__mutex_processes_waiting_queue.acquire()
        process_tuple2return = None, None, None
        if len(self.__waiting_list) > 0:
            process_tuple2return = self.__waiting_list.pop(0)
        self.__mutex_processes_waiting_queue.release()
        logger.info(f'process2return = {process_tuple2return}')
        return process_tuple2return
    
    def get_job_state(self, process_id: str, job_prefix: str):
        state2return = None
        if process_id in self.__processes_state_dict:
            state2return = self.__processes_state_dict[process_id].get_job_state(job_prefix)
        else:
            self.__mutex_processes_waiting_queue.acquire()
            for process_tuple in self.__waiting_list:
                if process_id in process_tuple:
                    state2return = State.Waiting
                    break
            self.__mutex_processes_waiting_queue.release()
        logger.info(f'state2return = {state2return} job_prefix = {job_prefix}')
        return state2return