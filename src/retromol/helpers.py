import errno 
import functools
import os 
import signal
import typing as ty

class TimeoutError(Exception):
    pass 

def timeout(seconds: int = 5, error_message: str = os.strerror(errno.ETIME)) -> ty.Callable:
    """
    Decorator to raise a TimeoutError when runtime exceeds the specified time.

    Parameters
    ----------
    seconds : int, optional
        Timeout in seconds, by default 5.
    error_message : str, optional
        Error message to be raised, by default os.strerror(errno.ETIME).
    
    Note: Timer only works on UNIX systems; also not thread-safe.
    """
    def decorator(func: ty.Callable) -> ty.Callable:

        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try: 
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result
        
        return wrapper
    
    return decorator