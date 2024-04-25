import multiprocessing

from .cli import app

if __name__ == "__main__":
    multiprocessing.freeze_support()
    app()
