from __future__ import print_function

import yaml
import argparse
import logging
import sevenbridges as sb
from sevenbridges.errors import SbgError

global config

try:
    fp = open('config.yaml')
    config = yaml.load(fp)
except:
    raise SbgError('config.yaml file missing!')

logger = logging.getLogger('star')
log_handler = logging.FileHandler(config['log_file'])
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
log_handler.setFormatter(formatter)
logger.addHandler(log_handler)
logger.setLevel(logging.DEBUG)

max_task_number = config['task_max_per_run']
project = config['project']
app = config['app']


def create_draft_tasks(api):
    logger.info('Creating draft tasks!')
    inputs = config['inputs']

    tar_list = list(api.files.query(project=project).all())
    tar_inputs = [_file for _file in tar_list if
                  _file.name.lower().endswith(inputs['input_archive_file'])]

    sjdbGTFfile_list = api.files.query(project=project,
                                       names=[inputs['sjdbGTFfile']])
    if len(sjdbGTFfile_list) == 0:
        raise SbgError(
            "Got 0 inputs! Expected at least 1 for input {}".format(
                config['sjdbGTFfile']
            ))

    genome_list = api.files.query(project=project,
                                           names=[inputs['genome']])
    if len(genome_list) == 0:
        raise SbgError(
            "Got 0 inputs! Expected at least 1 for input {}".format(
                config['genome']
            ))
    sjdbGTFfile = sjdbGTFfile_list[0]
    genomeFastaFile = genome_list[0]

    for _file in tar_inputs:
        task_name = 'Task.{}'.format(_file.name)
        inputs = {
            'sjdbGTFfile': [sjdbGTFfile],
            'genome': genomeFastaFile,
            'input_archive_file': _file
        }
        try:
            api.tasks.create(
                task_name, project, app, inputs=inputs, description=task_name
            )
        except SbgError as e:
            logger.error("Draft task was not created!", exc_info=e)
            raise SbgError("Draft task was not created!")


def run_tasks(api):
    logger.info('Running tasks!')

    running_tasks = list(
        api.tasks.query(project=project, limit=100, status='RUNNING').all()
    )
    queued_tasks = list(
        api.tasks.query(project=project, limit=100, status='QUEUED').all()
    )
    if len(running_tasks) + len(queued_tasks) >= max_task_number:
        logger.info("Maximum number of active tasks reached!")
        raise SbgError(
            'Unable to run! You already have {active} active tasks. '
            'Please try later!'.format
            (active=len(running_tasks) + len(queued_tasks)))

    draft_tasks = list(
        api.tasks.query(project=project, limit=100, status='DRAFT').all()
    )
    if len(draft_tasks) == 0:
        print('No draft tasks left to be run!')
        return

    executable_tasks = draft_tasks[0:max_task_number - len(running_tasks)]
    for task in executable_tasks:
        try:
            task.run()
        except SbgError as e:
            logger.error("Task was not started! Error happened ", exc_info=e)
            raise SbgError('Task was not started! Error happened')
        if task.status == 'DRAFT':
            logger.error("Task was not started! Task state is DRAFT!")
            raise SbgError("Task was not started! Task state is DRAFT!")


def status(api):
    logger.info('Fetching task statuses!')
    queued = api.tasks.query(project=project, status='QUEUED').total
    running = api.tasks.query(project=project, status='RUNNING').total
    completed = api.tasks.query(project=project, status='COMPLETED').total
    draft = api.tasks.query(project=project, status='DRAFT').total
    failed = api.tasks.query(project=project, status='FAILED').total
    aborted = api.tasks.query(project=project, status='ABORTED').total
    print("Draft={}, Queued={}, Running={}, Completed={},"
          " Failed={}, Aborted={} ".format(draft, queued,
                                           running, completed,
                                           failed, aborted)
          )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "option", nargs="?", help="create|run|status"
    )
    args = parser.parse_args()

    sb_config = sb.Config(url=config['api-url'], token=config['token'])
    api = sb.Api(config=sb_config)

    if args.option == 'create':
        create_draft_tasks(api)

    if args.option == 'run':
        run_tasks(api)

    if args.option == 'status':
        status(api)
