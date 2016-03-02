import logging

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > len(logging_levels):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)
