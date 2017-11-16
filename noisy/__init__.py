import os
os.environ['MPLBACKEND'] = 'Agg'

import pkg_resources
try:
    __version__ = pkg_resources.require("noisy")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"
