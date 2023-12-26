#!/usr/bin/env python
'''
setup.py file for GWforge package
'''
from setuptools import setup, find_packages
import os

def find_files(dirname, relpath=None, extensions=(".py", ".pyc"), ignore_dirs=None, ignore_files=None):
    def find_paths(directory):
        items = []
        for root, dirs, files in os.walk(directory):
            # Exclude specified directories if ignore_dirs is provided
            if ignore_dirs is not None:
                dirs[:] = [d for d in dirs if d not in ignore_dirs]

            for fname in files:
                path = os.path.join(root, fname)
                # Exclude specified files and extensions if ignore_files is provided
                if ignore_files is not None and (fname in ignore_files or any(path.endswith(ext) for ext in extensions)):
                    continue

                items.append(path)
        return items

    items = find_paths(dirname)
    relpath = relpath or dirname

    return [os.path.relpath(path, relpath) for path in items]

ignore_dirs = ['.git', '__pycache__/', 'profile_default/', '.ipynb_checkpoints', '.vscode', 'venv/', 'env/', 'virtualenv/']
ignore_files = ['.DS_Store', '*.pyc', 'Thumbs.db', '*.sqlite', '*.swp', '*.env', '.env']

setup(
    scripts=find_files('bin/', relpath='./', ignore_dirs=ignore_dirs, ignore_files=ignore_files),
    packages=find_packages()
)
