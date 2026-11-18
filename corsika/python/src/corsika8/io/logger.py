import logging

c8_logger = logging.getLogger("corsika")
if not c8_logger.hasHandlers():
    fmt = "[%(levelname)s] - %(name)s - %(message)s"
    myFormatter = logging.Formatter(fmt)
    handler = logging.StreamHandler()
    handler.setFormatter(myFormatter)
    c8_logger.addHandler(handler)
