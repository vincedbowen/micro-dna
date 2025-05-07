import functools
import threading
import time

def status_message(message: str, delay: float = 0.5):
    """Decorator to print a status message with a loading animation.

    Args:
        message (str): Message to display while the function is running.
        delay (float, optional): delay for dot animation. Defaults to 0.5.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            stop_event = threading.Event()
            def print_dots():
                while not stop_event.is_set():
                    print('.', end='', flush=True)
                    time.sleep(delay)

            print(f"{message} ", end='', flush=True)
            dot_thread = threading.Thread(target=print_dots)
            dot_thread.start()

            try:
                result = func(*args, **kwargs)
            finally:
                stop_event.set()
                dot_thread.join()
                print("✅")
            return result

        return wrapper
    return decorator

