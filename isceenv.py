#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper to run ISCE-environment commands (topsApp.py, convert_isce, gdal2isce)
from a process running in the base conda environment. The Python workflow runs
in base; only the isce-dependent shell commands are dispatched into the
isce2grimp env via `conda run`, so the user never has to activate it.

Commands are teed to a run log (if set) and raise on non-zero exit so a failed
ISCE step aborts the pipeline instead of silently continuing.
"""
import os
import shutil
from subprocess import call

_ISCE_ENV = 'isce2grimp'
_LOG = None


class IsceCommandError(RuntimeError):
    ''' An ISCE command exited non-zero. Subclasses RuntimeError so existing
    handlers still catch it, but carries the exit code so the caller can
    report it. '''

    def __init__(self, command, returncode):
        super().__init__(f'isceShell: command failed '
                         f'(exit {returncode}): {command}')
        self.command = command
        self.returncode = returncode


def condaExe():
    ''' Absolute path to conda so the shell-out does not depend on PATH. '''
    return (os.environ.get('CONDA_EXE') or shutil.which('conda')
            or '/home/ian/miniforge3/bin/conda')


def setEnv(name):
    ''' Override the conda env used for ISCE commands (from project yaml). '''
    global _ISCE_ENV
    if name:
        _ISCE_ENV = name


def getEnv():
    return _ISCE_ENV


def setLog(path):
    ''' Tee ISCE command output into this log file (None disables). '''
    global _LOG
    _LOG = path


def isceShell(command, ompThreads=None, ompPlaces=None, check=True):
    ''' Run a command inside the isce2grimp conda env. OMP settings are
    exported first so topsApp picks them up. Output is teed to the run log;
    a non-zero exit raises IsceCommandError (a RuntimeError carrying the exit
    code) unless check=False.

    Uses bash with pipefail so the tee pipe does not mask the command's exit
    code. '''
    prefix = ''
    if ompThreads is not None:
        prefix += f'OMP_NUM_THREADS={ompThreads} '
    if ompPlaces is not None:
        prefix += f"OMP_PLACES='{ompPlaces}' "
    full = (f'set -o pipefail; {prefix}{condaExe()} run -n {_ISCE_ENV} '
            f'--no-capture-output {command}')
    if _LOG is not None:
        full += f' 2>&1 | tee -a {_LOG}'
    print(full)
    rc = call(full, shell=True, executable='/bin/bash')
    if check and rc != 0:
        raise IsceCommandError(command, rc)
    return rc
