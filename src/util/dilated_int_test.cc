/// @brief  Z-order or Morton style indexing utilities for arbitrary dimension
/// @author will sheffler

// inspired by code from here: http://www.marcusbannerman.co.uk/dynamo
// see "Converting to and from Dilated Integers"(doi:10.1109/TC.2007.70814)

#include <util/dilated_int.hh>
#include <gtest/gtest.h>

namespace scheme {
namespace util {

template<uint64_t D> bool test_dilated_int() {
	uint64_t maxval = util::impl::MAX_DILATABLE<D>::VAL;
	maxval = std::min(maxval, uint64_t(2097151));
	for (uint64_t i(0); i < maxval; ++i){
		uint64_t dilated = util::dilate<D>(i);
		uint64_t undilated = util::undilate<D>(dilated);
		EXPECT_EQ(undilated,i);
    }
    return true;
}

TEST(DILATED_INT,test_64bit_DIM01){ util::test_dilated_int< 1>(); }
TEST(DILATED_INT,test_64bit_DIM02){ util::test_dilated_int< 2>(); }
TEST(DILATED_INT,test_64bit_DIM03){ util::test_dilated_int< 3>(); }
TEST(DILATED_INT,test_64bit_DIM04){ util::test_dilated_int< 4>(); }
TEST(DILATED_INT,test_64bit_DIM05){ util::test_dilated_int< 5>(); }
TEST(DILATED_INT,test_64bit_DIM06){ util::test_dilated_int< 6>(); }
TEST(DILATED_INT,test_64bit_DIM07){ util::test_dilated_int< 7>(); }
TEST(DILATED_INT,test_64bit_DIM08){ util::test_dilated_int< 8>(); }
TEST(DILATED_INT,test_64bit_DIM09){ util::test_dilated_int< 9>(); }
TEST(DILATED_INT,test_64bit_DIM10){ util::test_dilated_int<10>(); }
TEST(DILATED_INT,test_64bit_DIM11){ util::test_dilated_int<11>(); }
TEST(DILATED_INT,test_64bit_DIM12){ util::test_dilated_int<12>(); }
TEST(DILATED_INT,test_64bit_DIM13){ util::test_dilated_int<13>(); }
TEST(DILATED_INT,test_64bit_DIM14){ util::test_dilated_int<14>(); }
TEST(DILATED_INT,test_64bit_DIM15){ util::test_dilated_int<15>(); }
TEST(DILATED_INT,test_64bit_DIM16){ util::test_dilated_int<16>(); }
TEST(DILATED_INT,test_64bit_DIM17){ util::test_dilated_int<17>(); }
TEST(DILATED_INT,test_64bit_DIM18){ util::test_dilated_int<18>(); }
TEST(DILATED_INT,test_64bit_DIM19){ util::test_dilated_int<19>(); }
TEST(DILATED_INT,test_64bit_DIM20){ util::test_dilated_int<20>(); }
TEST(DILATED_INT,test_64bit_DIM21){ util::test_dilated_int<21>(); }
TEST(DILATED_INT,test_64bit_DIM22){ util::test_dilated_int<22>(); }
TEST(DILATED_INT,test_64bit_DIM23){ util::test_dilated_int<23>(); }
TEST(DILATED_INT,test_64bit_DIM24){ util::test_dilated_int<24>(); }
TEST(DILATED_INT,test_64bit_DIM25){ util::test_dilated_int<25>(); }
TEST(DILATED_INT,test_64bit_DIM26){ util::test_dilated_int<26>(); }
TEST(DILATED_INT,test_64bit_DIM27){ util::test_dilated_int<27>(); }
TEST(DILATED_INT,test_64bit_DIM28){ util::test_dilated_int<28>(); }
TEST(DILATED_INT,test_64bit_DIM29){ util::test_dilated_int<29>(); }
TEST(DILATED_INT,test_64bit_DIM30){ util::test_dilated_int<30>(); }
TEST(DILATED_INT,test_64bit_DIM31){ util::test_dilated_int<31>(); }
TEST(DILATED_INT,test_64bit_DIM32){ util::test_dilated_int<32>(); }
// above 32 makes no sense for 64 bit



}
}
