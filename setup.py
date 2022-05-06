from setuptools import setup
import re


with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

with open('pepmatch/version.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

setup(
    name='pepmatch',
    version=version,
    description='Peptide and epitope search against a reference proteome with specified mismatches.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/IEDB/PEPMatch',
    author='Daniel Marrama',
    author_email='dmarrama@lji.org',
    packages=['pepmatch'],
    install_requires=['numpy>=1.18',
                      'pandas>=1.1',
                      'biopython>=1.5',
                      'python-Levenshtein>=0.11',
                      'openpyxl>=3.0.0'],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pepmatch-preprocess = pepmatch.preprocessor:run',
            'pepmatch-match = pepmatch.matcher:run'
        ]
    }
)
