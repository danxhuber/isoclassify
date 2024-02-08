import os
import importlib.metadata

try:
     __version__ = importlib.metadata.version(__package__)  # only works if package installed via pip
except:
     print('')

# Absolute path to the package directory
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Path to external isoclassify data (e.g. MESA grid)
# TODO: Put all this in a utility function
DATADIR = None

if 'ISOCLASSIFY' in os.environ:
    DATADIR = os.environ['ISOCLASSIFY']
else:
    DATADIR = os.path.join(os.path.expanduser('~'), '.isoclassify')

if not os.path.isdir(DATADIR):
    # If DATADIR doesn't exist, make a new directory to accommodate data
    try:
        os.mkdir(DATADIR)
    except OSError:
        # If cannot make new directory, we assume it to be CWD
        DATADIR = '.'
