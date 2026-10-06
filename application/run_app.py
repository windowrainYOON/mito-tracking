"""Launcher used by PyInstaller and for running from source: python run_app.py"""
import multiprocessing
import sys

from mito_app.app import main

if __name__ == '__main__':
    multiprocessing.freeze_support()
    sys.exit(main())
