from setuptools import setup, Extension

hamming_module = Extension(
    'pepmatch.hamming',
    sources=['pepmatch/hamming.c']
)

setup(
    ext_modules=[hamming_module]
)
