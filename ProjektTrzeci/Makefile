docs:
	doxygen DoxygenConfig

build-sync:
	f2py -c -m Matrix SyncMatrix.f90

run-python:
	python Pyapp.py

build-async:
	ifort -coarray AsyncMatrix.f90 -o AMatrix

run-async:
	./AMatrix

time-compare:
	time python Pyapp.py
	time ./AMatrix