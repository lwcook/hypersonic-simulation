#!/bin/bash

rm -r build

nosetests --with-coverage --cover-package=horsetailmatching --cover-html --cover-erase

coverage-badge -o tests/coverage.svg

python setup.py sdist upload

git add -A
git commit -m "Publishing latest version"
git push
