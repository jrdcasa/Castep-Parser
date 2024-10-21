import logging
from datetime import datetime
import sys
from setuptools import setup


# Formatter for the logger
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# Main setup
if __name__ == '__main__':

    # Creating the logger to install.log file ==================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Starting installation!!!! at {}\n\n".format(nowm)
    m1 += "\t\tPython version {}.{}\n".format(sys.version_info.major, sys.version_info.minor)
    print(m1) if logger is None else logger.info(m1)

    m1 = "\t\t ============== SYS PATH ==============\n"
    for item in sys.path:
        m1 += "\t\t\t" + item + "\n"
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 += "\t\t ============== END SYS PATH ==============\n"
    print(m1) if logger is None else logger.info(m1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t ============== RUNNING SETUP FROM SETUPTOOLS {} ==============\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    setup()

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Installation Done!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)