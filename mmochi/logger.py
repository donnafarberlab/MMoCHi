import logging
import sys

logg = logging.getLogger(__name__)
level_names = dict(INFO=logging.INFO, DEBUG=logging.DEBUG, WARNING=logging.WARNING, ERROR=logging.ERROR)

def initiate_log(file_name, file_level='DEBUG', stream_level='INFO', file_mode='w'):
    '''Creates a log files with the name log_file, with levels defined by file level and stream level.
    levels can be one of: "DEBUG", "INFO", "WARNING","ERROR" '''
    _file_level = level_names[file_level]
    _stream_level = level_names[stream_level]
    streamHandler = logging.StreamHandler(sys.stderr)
    streamHandler.setFormatter(logging.Formatter("%(message)s"))#("%(levelname)s - %(message)s"))
    streamHandler.setLevel(logging.INFO)
    logging.basicConfig(format="%(levelname)s - %(message)s",level=_stream_level)
    logg.propagate = False
    if (logg.hasHandlers()):
        logg.handlers.clear()
    if not file_name is None:
        fileHandler = logging.FileHandler('Classifier.log',mode=file_mode)
        fileHandler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        fileHandler.setLevel(_file_level)
        logg.addHandler(fileHandler)
    logg.addHandler(streamHandler)
    logg.setLevel(logging.DEBUG)
    if not file_name is None:
        logg.debug(f'Set up logger with file_level {file_level} and stream_level {stream_level}')
    else:
        logg.debug(f'Set up logger with stream_level {stream_level}')
    return

