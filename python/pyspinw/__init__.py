from __future__ import annotations

__author__ = "github.com/wardsimon"
__version__ = "0.0.0"

import os
import libpymcr
from . import plotting


# Generate a list of all the MATLAB versions available
_VERSION_DIR = os.path.join(os.path.dirname(__file__), 'ctfs')
_VERSIONS = []
for file in os.scandir(_VERSION_DIR):
    if file.is_file() and file.name.endswith('.ctf'):
        _VERSIONS.append({'file':    os.path.join(_VERSION_DIR, file.name),
                          'version': 'R' + file.name.split('.')[0].split('SpinW_')[1]
                          })

VERSION = ''
INITIALIZED = False


class Matlab(libpymcr.Matlab):
    def __init__(self, matlab_path: Optional[str] = None, matlab_version: Optional[str] = None):
        """
        Create a MATLAB instance with the correct compiled library for the MATLAB version specified. If no version is
        specified, the first version found will be used. If no MATLAB versions are found, a RuntimeError will be
        raised. If a version is specified, but not found, a RuntimeError will be raised.

        :param matlab_path: Path to the root directory of the MATLAB installation or MCR installation.
        :param matlab_version: Used to specify the version of MATLAB if the matlab_path is given or if there is more
        than 1 MATLAB installation.
        """

        global INITIALIZED
        global VERSION
        if INITIALIZED:
            super().__init__(VERSION, mlPath=matlab_path)
        elif matlab_version is None:
            for version in _VERSIONS:
                try:
                    print(f"Trying MATLAB version: {version['version']} ({version['file']}))")
                    super().__init__(version['file'], mlPath=matlab_path)
                    INITIALIZED = True
                    VERSION = version['version']
                except RuntimeError:
                    continue
        else:
            ctf = [version for version in _VERSIONS if version['version'].lower() == matlab_version.lower()]
            if len(ctf) == 0:
                raise RuntimeError(
                    f"Compiled library for MATLAB version {matlab_version} not found. Please use: [{', '.join([version['version'] for version in _VERSIONS])}]\n ")
            else:
                ctf = ctf[0]
            try:
                super().__init__(ctf['file'], mlPath=matlab_path)
                INITIALIZED = True
                VERSION = ctf['version']
            except RuntimeError:
                pass
        if not INITIALIZED:
            raise RuntimeError(
                f"No supported MATLAB versions [{', '.join([version['version'] for version in _VERSIONS])}] found.\n "
                f"If installed, please specify the root directory (`matlab_path` and `matlab_version`) of the MATLAB "
                f"installation.")

from .matlab_wrapped import *
