# distutils: language = c++

cdef extern from "dilated_int.hh" namespace "util":
	cdef bint test_dilated_int()


def test_zorder():
	return test_dilated_int()
