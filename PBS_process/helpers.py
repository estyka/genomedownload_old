import functools
from datetime import datetime

def timeit(func):
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        start_time = datetime.now()
        result = func(*args, **kwargs)
        elapsed_time = datetime.now() - start_time
        print('function [{}] finished in {} ms'.format(
            func.__name__, elapsed_time))
        return result

    return new_func