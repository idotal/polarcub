all:
	python3 setup.py build_ext --inplace

clean:
	rm -rf *.c *.so *.html VectorDistributions/*.c  VectorDistributions/*.so VectorDistributions/*.html build
