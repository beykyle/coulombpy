build:
	python3 -m numpy.f2py -c coulfg.pyf coulfg1.f 


conf:
	python3 -m numpy.f2py coulfg1.f -m coulfg -h coulfg.pyf

clean:
	rm *.so 

confclean:
	rm *.pyf
