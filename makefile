conf:
	python3 -m numpy.f2py coulfg1.f -m coulfg -h coulfg.pyf

build:
	python3 -m numpy.f2py -c coulfg.pyf coulfg1.f 

clean:
	rm *.so 

confclean:
	rm *.pyf
