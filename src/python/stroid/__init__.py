import io
import sys

from ._stroid import *

from ._stroid import config
from ._stroid import exceptions
from ._stroid import IO
from ._stroid import refinement
from ._stroid import stats
from ._stroid import GenerateMesh
from ._stroid import StroidMesh

sys.modules['stroid.config'] = config
sys.modules['stroid.exceptions'] = exceptions
sys.modules['stroid.IO'] = IO
sys.modules["stroid.refinement"] = refinement
sys.modules["stroid.stats"] = stats


__all__ = ['config', 'exceptions', 'IO', 'refinement', 'stats', 'GenerateMesh', 'StroidMesh']

import importlib.metadata


try:
    _meta = importlib.metadata.metadata('stroid')
    __version__ = _meta['Version']
    __license__ = _meta['License']
    __description__ = _meta['Summary']
    __author__ = 'Emily M. Boudreaux'
    __url__ = 'https://github.com/4D-STAR/stroid'
except importlib.metadata.PackageNotFoundError :
    __version__ = 'unknown - Package not installed'
    __license__ = 'GNU General Public License v3.0'
    __email__ = 'emily.boudreaux@dartmouth.edu'
    __url__ = 'https://github.com/4D-STAR/stroid'


import os
from pathlib import Path
from typing import List

_PACKAGE_DIR = Path(__file__).resolve().parent

