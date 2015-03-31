gnu_test: 
	g++ gnu_test.cc -o gnu_test -lm -lgsl -lgslcblas -std=c++11 -I/soft/gsl/gnu/1.15/include/ -L/soft/gsl/gnu/1.15/lib/ -fopenmp

serial_search: 
	CC serial_search.cc error_function.cc -o serial_search -lm -lgsl -lgslcblas -std=c++11 -I/soft/gsl/gnu/1.15/include/ -L/soft/gsl/gnu/1.15/lib/ -fopenmp
