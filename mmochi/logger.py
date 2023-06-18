import logging
import sys

def _addLoggingLevel(levelName, levelNum, methodName=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `levelName` becomes an attribute of the `logging` module with the value
    `levelNum`. `methodName` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `methodName` is not specified, `levelName.lower()` is
    used.

    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present 

    Example
    -------
    >>> _addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not methodName:
        methodName = levelName.lower()

    if hasattr(logging, levelName):
        raise AttributeError('{} already defined in logging module'.format(levelName))
    if hasattr(logging, methodName):
        raise AttributeError('{} already defined in logging module'.format(methodName))
    if hasattr(logging.getLoggerClass(), methodName):
        raise AttributeError('{} already defined in logger class'.format(methodName))

    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(levelNum):
            self._log(levelNum, message, args, **kwargs)
    def logToRoot(message, *args, **kwargs):
        logging.log(levelNum, message, *args, **kwargs)

    logging.addLevelName(levelNum, levelName)
    setattr(logging, levelName, levelNum)
    setattr(logging.getLoggerClass(), methodName, logForLevel)
    setattr(logging, methodName, logToRoot)

_addLoggingLevel('PRINT',logging.DEBUG+5)

logg = logging.getLogger(__name__)
level_names = dict(INFO=logging.INFO, PRINT=logging.PRINT, DEBUG=logging.DEBUG, 
                   WARNING=logging.WARNING, ERROR=logging.ERROR)

class _OneLevelStreamHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno == self.level:
            super().emit(record)

class _OneLevelFileHandler(logging.FileHandler):
    def emit(self, record):
        if record.levelno == self.level:
            super().emit(record)

def _initiate_log():
    '''Creates a log files with the name log_file'''
    logging.basicConfig(format="%(levelname)s - %(message)s",level=logging.INFO)
    logg.propagate = False
    if (logg.hasHandlers()):
        logg.handlers.clear()
    errStreamHandler = logging.StreamHandler(sys.stderr)
    errStreamHandler.setFormatter(logging.Formatter("%(message)s"))
    errStreamHandler.setLevel(logging.INFO)
    logg.addHandler(errStreamHandler)
    
    stoutStreamHandler = _OneLevelStreamHandler(sys.stdout)
    stoutStreamHandler.setFormatter(logging.Formatter("%(message)s"))
    stoutStreamHandler.setLevel(logging.PRINT)
    logg.addHandler(stoutStreamHandler)
    
    logg.setLevel(logging.DEBUG)
    logg.debug('Set up logger.')
    return

def log_to_file(file_name, file_mode='w'):
    '''
    Enable logging for all functions in the MMoCHi package
    
    Parameters
    ----------
    file_name
        The file path and file name (without a ".log" suffix) detailing where a log file should be saved.
    file_mode
        Either "w" to overwrite a preexisting file, or "a" to append logs to a preexisting file. 
    '''
    fileHandler = logging.FileHandler(f'{file_name}.log',mode=file_mode)
    fileHandler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    fileHandler.setLevel(logging.DEBUG)
    logg.addHandler(fileHandler)