from setuptools import setup, Extension
import re


with open('README.md', 'r', encoding='utf-8') as f:
  long_description = f.read()


with open('pepmatch/version.py', 'r') as f:
  version = re.search(
    r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    f.read(),
    re.MULTILINE).group(1)

hamming_module = Extension(
  'pepmatch.hamming',
  sources=['pepmatch/hamming.c']
)

setup(
  name='pepmatch',
  version=version,
  description='Search tool for peptides and epitopes within a proteome, '
              'while considering potential residue substitutions.',
  long_description=long_description,
  long_description_content_type='text/markdown',
  url='https://github.com/IEDB/PEPMatch',
  author='Daniel Marrama',
  author_email='dmarrama@lji.org',
  packages=['pepmatch'],
  ext_modules=[hamming_module],
  install_requires=[
    'numpy>=1.18',
    'pandas>=1.1',
    'biopython>=1.5',
    'python-Levenshtein>=0.11',
    'openpyxl>=3.0.0'
  ],
  zip_safe=False,
  entry_points={
    'console_scripts': [
      'pepmatch-preprocess = pepmatch.shell:run_preprocessor',
      'pepmatch-match = pepmatch.shell:run_matcher'
    ]
  }
)
