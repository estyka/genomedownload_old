from flask import Flask, render_template, request, url_for, redirect
import uuid
from SharedConsts import UI_CONSTS
from Job_Manager_API import Job_Manager_API
import os

# from ..utils import logger
print(__name__)
app = Flask(__name__)

MAX_NUMBER_PROCESS = 3
UPLOAD_FOLDERS_ROOT_PATH = '/bioseq/data/results/genomedownload/'

process_id2update = []

def update_html(process_id, state):
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
