import multiprocessing

from .anglerfish import anglerfish

if __name__ == "__main__":
    multiprocessing.freeze_support()
    anglerfish()
