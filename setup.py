from setuptools import setup

# TODO: option to install BC and/or MIST grids on install

# Load version
__version__ = None
# exec(open('interstellar/version.py').read())

# Load requirements
requirements = None
# with open('requirements.txt') as file:
#     requirements = file.read().splitlines()

setup(
    name='isoclassify',
    version=__version__,
    description='Isoclassify',
    packages=['isoclassify'],
    author='Daniel Huber',
    install_requires=requirements,
)
