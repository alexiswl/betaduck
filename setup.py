from setuptools import setup, find_packages
import sys

import betaduck.version

if sys.version_info < (3, 7):
    sys.exit('Sorry, Python < 3.7 is not supported')

with open("requirements.txt", 'r') as file_h:
    requirements = [line.strip() for line in file_h.readlines()]

long_description = """

Betaduck is a set of tools for handling, basecalling and producing quality metrics of Promethion beta data.

"""

setup(
    name='betaduck',
    version=betaduck.__version__,
    packages=find_packages(),
    provides=['betaduck'],
    requires=['python (>=3.7)'],
    install_requires=requirements,
    url='github.com/alexiswl/poreduck',
    license='GPL',
    author='Alexis Lucattini',
    author_email='alexis.lucattini@agrf.org.au',
    description='Nanopore data handling for Promethion Beta.',
    package_dir={'betaduck': "betaduck"},
    entry_points={
        'console_scripts': ['betaduck=betaduck.betaduck_main:main']
    }
)
