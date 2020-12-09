from setuptools import setup


with open('README.md', 'r') as f:
    long_description = f.read()


setup(
    name='pepmatch',
    version='0.3',
    description='Peptide/epitope search against a reference proteome with specified mismatches.',
    long_description=long_description,
    url='https://gitlab.lji.org/dmarrama/pepmatch',
    author='Daniel Marrama',
    author_email='dmarrama@lji.org',
    packages=['pepmatch'],
    install_requires=['numpy>=1.18',
                  'pandas>=1.1',
                  'biopython>=1.5',
                  'python-Levenshtein>=0.11'],
    zip_safe=False
)
