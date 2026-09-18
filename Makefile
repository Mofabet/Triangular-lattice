.PHONY: install test fast lint figures bench clean demo

install:
	pip install -e ".[dev,fast]"

test:
	pytest -q

fast:                       ## skip the long statistical-sampling tests
	pytest -q -m "not slow" -x

figures:                    ## regenerate every figure (~9 min on one core)
	python examples/00_potential.py
	python examples/01_validation.py
	python examples/02_melting.py heat
	python examples/02_melting.py cool
	python examples/03_figures.py
	python examples/04_defects.py
	python examples/05_lattices.py

bench:
	trilattice bench

demo:                       ## the original input deck, done properly
	trilattice run legacy/start.txt -T 300 -n 20000 --xyz demo.xyz

clean:
	rm -rf build dist .pytest_cache **/__pycache__ *.egg-info
