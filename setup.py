from setuptools import setup

# Load version
__version__ = None
# exec(open('isoclassify/version.py').read())

# Load requirements
requirements = None
with open('requirements.txt') as file:
    requirements = file.read().splitlines()

# mwdust requires manual install to download maps
requirements.remove('mwdust')

setup(
    name='isoclassify',
    version=__version__,
    description='Isoclassify',
    packages=['isoclassify'],
    author='Daniel Huber',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'isoclassify = isoclassify.isoclassify:main',
        ],
    }
)
