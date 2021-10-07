.PHONY: all
all:
	python3 setup.py build_ext --inplace


.PHONY: test
test:
	pytest -sv .


.PHONY: clean
clean:
	rm -rf build permanent/*.so
