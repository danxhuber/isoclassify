import os
from setuptools import setup

# Set version (increment this when making a new release)
version = "1.4"

# Load requirements
requirements = None
with open('requirements.txt') as file:
    requirements = file.read().splitlines()

# mwdust requires manual install to download maps
# TODO: add warning if `mwdust` is not installed
requirements.remove('mwdust')

# Configure data directory
if 'ISOCLASSIFY' in os.environ:
    # For those who've already configured environment variables
    datadir = os.environ['ISOCLASSIFY']
else:
    # Create isoclassify directory in user home
    datadir = os.path.join(os.path.expanduser('~'), '.isoclassify')

if not os.path.isdir(datadir):
    # If datadir doesn't exist, make a new one upon install
    os.mkdir(datadir)

# Package description
desc = 'Python codes to perform stellar classifications given any set of input observables.'

setup(
    name='isoclassify',
    version=version,
    description=desc,
    package_dir={
        'isoclassify': 'isoclassify', 
        'isoclassify.direct': 'isoclassify/direct', 
        'isoclassify.grid': 'isoclassify/grid',
        },
    packages=['isoclassify', 'isoclassify.direct', 'isoclassify.grid'],
    include_package_data=True,  # <-- includes 'isoclassify/data'
    author='Daniel Huber',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'isoclassify = isoclassify.isoclassify:main',
        ],
    }
)
