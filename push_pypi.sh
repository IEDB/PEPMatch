#!/bin/bash

rm -rf dist
ipython3 setup.py sdist bdist_wheel
twine upload dist/*
