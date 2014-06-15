# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector

#cdef extern from "Nest.hh" namespace "Eigen":
#    cdef cppclass Array[Float,DIM,1]

cdef extern from "NEST.hh" namespace "scheme::nest":
    cdef void test()


test()

