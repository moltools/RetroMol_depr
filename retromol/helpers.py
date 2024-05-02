"""This module contains helper functions and decorators for the retromol
package.
"""

import errno
import functools
import os
import signal
import typing as ty


class RetroMolTimeout(Exception):
    """Custom exception to raise when a function exceeds the specified time."""

    pass


def timeout(
    seconds: int = 5, error_message: str = os.strerror(errno.ETIME)
) -> ty.Callable:
    """Raise a RetroMolTimeout when runtime exceeds the specified time.

    :param seconds: The number of seconds before the function times out.
    :type seconds: int
    :param error_message: The error message to display when the function times
        out.
    :type error_message: str
    :return: Decorator.
    :rtype: ty.Callable
    """

    def decorator(func: ty.Callable) -> ty.Callable:

        def _handle_timeout(signum, frame):
            raise RetroMolTimeout(error_message)

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
