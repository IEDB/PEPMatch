from setuptools import setup


with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()


setup(
    name='pepmatch',
    version='0.7.15',
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
                  'python-Levenshtein>=0.11'],
    zip_safe=False
)
