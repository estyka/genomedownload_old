from flask import Flask, render_template, request, url_for, redirect, Response, jsonify, send_file
import uuid
from SharedConsts import UI_CONSTS
from Job_Manager_API import Job_Manager_API
import os
from utils import State, logger
import time


# from ..utils import logger
print(__name__)
app = Flask(__name__)

MAX_NUMBER_PROCESS = 3
UPLOAD_FOLDERS_ROOT_PATH = '/bioseq/data/results/genomedownload/'
TIME_OF_STREAMING_UPDATE_REQUEST_BEFORE_DELETING_IT_SEC = 1200


process_id2update = []

def update_html(process_id, state):
    #TODO: see that the states are correct in every stage (running, waiting, finished...)
    logger.info(f'process_id = {process_id} state = {state}')
    if process_id:
        process_id2update.append(process_id)

@app.route('/remove_update/<process_id>')
def remove_update(process_id):
    logger.info(f'process_id = {process_id}')
    if process_id in process_id2update:
        logger.info(f'removing {process_id} from process_id2update')
        process_id2update.remove(process_id)
    return jsonify('data')

@app.route('/stream/')
def stream():
    # function to stream data to client
    requests_time_dict = {}
    TIME_BETWEEN_BROADCASTING_EVENTS = 0.1
    
    def eventStream():
        while True:
            if len(process_id2update):
                for process_id in process_id2update:
                    yield 'data: {}\n\n'.format(process_id)
                    time.sleep(TIME_BETWEEN_BROADCASTING_EVENTS)
                    time_broadcasting_process_event = requests_time_dict.get(process_id, 0)
                    time_broadcasting_process_event += TIME_BETWEEN_BROADCASTING_EVENTS
                    requests_time_dict[process_id] = time_broadcasting_process_event

                max_broadcasting_event = max(requests_time_dict, key=requests_time_dict.get)
                if requests_time_dict[max_broadcasting_event] >= TIME_OF_STREAMING_UPDATE_REQUEST_BEFORE_DELETING_IT_SEC:
                    logger.info(f'removing max_broadcasting_event = {process_id} as it reached the max amount of time broadcasting')
                    requests_time_dict.pop(max_broadcasting_event)
            else:
                requests_time_dict.clear()
                time.sleep(1)
    return Response(eventStream(), mimetype="text/event-stream")

manager = Job_Manager_API(MAX_NUMBER_PROCESS, UPLOAD_FOLDERS_ROOT_PATH, update_html)

@app.route('/process_state/<process_id>')
def process_state(process_id):
    job_state = manager.get_kraken_job_state(process_id)
    if job_state == None:
        return redirect(url_for('error', error_type=UI_CONSTS.UI_Errors.UNKNOWN_PROCESS_ID.name))
    if job_state != State.Finished:
        kwargs = {
            "process_id": process_id,
            "text": UI_CONSTS.states_text_dict[job_state],
            "message_to_user": UI_CONSTS.PROCESS_INFO_KR,
            "gif": UI_CONSTS.states_gifs_dict[job_state],
        }
        return render_template('process_running.html', **kwargs)
    else:
        return redirect(url_for('download_file', process_id=process_id))

@app.route('/download_file/<process_id>', methods=['GET', 'POST'])
def download_file(process_id):
    logger.info(f'request.method: {request.method}')
    if request.method == 'POST':
        dir2send = manager.export_dir(process_id)
        if dir2send == None:
            logger.warning(f'failed to export dir, process_id = {process_id}, dir2send = {dir2send}')
            return redirect(url_for('error', error_type=UI_CONSTS.UI_Errors.EXPORT_FILE_UNAVAILABLE.name)) #TODO: change error to be export_dir_unavailable
        logger.info(f'exporting, process_id = {process_id}, dir2send = {dir2send}')
        return send_file(dir2send, mimetype='application/octet-stream')
    return render_template('export_file.html')

@app.route('/error/<error_type>')
def error(error_type):
    # checking if error_type exists in error enum
    try:
        return render_template('error_page.html', error_text=UI_CONSTS.UI_Errors[error_type].value)
    except:
        return render_template('error_page.html', error_text=f'Unknown error, \"{error_type}\" is not a valid error code')

@app.route('/display_error/<error_text>')
def display_error(error_text):
    return render_template('error_page.html', error_text=error_text)

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        organism_name = request.form.get('text', None)
        email_address = request.form.get('email', None)
        if organism_name is not None and email_address is not None:
            new_process_id = manager.get_new_process_id()
            folder2save_file = os.path.join(UPLOAD_FOLDERS_ROOT_PATH, new_process_id)
            os.mkdir(folder2save_file)
            man_results = manager.add_process(new_process_id, email_address, organism_name)
            if not man_results:
                logger.warning(f'job_manager_api can\'t add process')
                return redirect(url_for('display_error', error_text=UI_CONSTS.ALERT_USER_TEXT_CANT_ADD_PROCESS))
            logger.info(f'process added man_results = {man_results}, redirecting url')
            return redirect(url_for('process_state', process_id=new_process_id))

        # if organism_name is None:
        #     #logger.warning(f'organism_name not available')
        #     return redirect(url_for('display_error', error_text=UI_CONSTS.ALERT_USER_TEXT_INVALID_MAIL)) #change
        # if email_address is None:
        #     #logger.warning(f'email_address not available')
        #     return redirect(url_for('display_error', error_text=UI_CONSTS.ALERT_USER_TEXT_INVALID_MAIL))
        # #logger.info(f'organism name = {organism_name}, email_address = {email_address}')
    return render_template('home.html')
    # return render_template("export_file.html")
    #return render_template("process_running.html")
