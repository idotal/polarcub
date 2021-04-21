all:
	rm -f __init__.py
	python3 setup.py build_ext --inplace
	touch __init__.py

clean:
	rm -f *.c *.so *.html __init__.py
	rm -rf build __pycache__
