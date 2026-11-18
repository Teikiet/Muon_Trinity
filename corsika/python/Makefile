##
# ##
# corsika
#
# @file

# find python3
PYTHON=`which python3`

# our testing targets
.PHONY: tests flake black isort all

all: black flake isort mypy tests

tests:
	${PYTHON} -m pytest --cov=src/corsika8 tests

flake:
	${PYTHON} -m flake8 src/corsika8 tests examples

black:
	${PYTHON} -m black -t py37 src/corsika8 tests examples

isort:
	${PYTHON} -m isort --atomic src/corsika8 tests examples

mypy:
	${PYTHON} -m mypy src/corsika8

# end
