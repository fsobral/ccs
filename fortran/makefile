gs_comp:
	gfortran gs_comp.f90 gentests.f90 gs.f90 -fopenmp -Wunused -O3 \
	-o ../bin/gs_comp

unittests:
	gfortran gentests.f90 unittests.f90 -llapack -o ../bin/unittests

clean:
	rm ../bin/* -fr
