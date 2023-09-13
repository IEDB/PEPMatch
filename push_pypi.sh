#!/bin/bash

rm -rf dist
python setup.py sdist bdist_wheel
auditwheel repair dist/pepmatch-0.9.5-cp311-cp311-linux_x86_64.whl
twine upload wheelhouse/*
