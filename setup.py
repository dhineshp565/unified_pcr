"""
Setup script for ONT Target Sequencing Pipeline
"""

from setuptools import setup, find_packages
import os


def read_file(filename):
    """Read file contents"""
    filepath = os.path.join(os.path.dirname(__file__), filename)
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return f.read()
    return ''


setup(
    name='ont_targseq',
    version='1.0.0',
    description='Unified pipeline for Oxford Nanopore Technologies targeted amplicon sequencing',
    long_description=read_file('docs/README.md'),
    long_description_content_type='text/markdown',
    author='unified_pcr',
    author_email='',
    url='https://github.com/dhineshp565/unified_pcr',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.7',
    install_requires=[
        # Core dependencies (minimal for basic functionality)
    ],
    extras_require={
        'dev': [
            'pytest>=7.0.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'ont-targseq=ont_targseq.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='ont nanopore sequencing bioinformatics amplicon targeted-sequencing',
)
