#include <gtest/gtest.h>

#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>

#include "scheme/util/SimpleArray.hh"

#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <map>

namespace scheme { namespace util {

using std::cout;
using std::endl;

template<class K, class V, class GetMap, size_t SEG=0>
struct SegmentedMap {
	typedef typename GetMap::template apply< K, util::SimpleArray<1<<SEG,V,true> >::type MAP;
	typedef typename MAP::const_iterator const_iterator;

	MAP map;

	SegmentedMap() { map.set_empty_key(std::numeric_limits<K>::max()); }

	const_iterator find(K const & k) const {
		return map.find(k>>SEG);
	}

	const_iterator end() const {
		return map.end();
	}

	void insert(std::pair<K,V> const & v){
		map.insert(std::make_pair(v.first>>SEG,v.second));
	}

	V const & operator[](K const & k) const {
		return map.find(k>>SEG)->second.operator[]( k % (1<<SEG) );
	}

	V & operator[](K const & k) {
		return map[k>>SEG][ k % (1<<SEG) ];
	}
};

template<class Map>
void test_map(Map & h){
	int64_t MAXIDX = 430l*1000*1000l;
	int64_t NFILL = MAXIDX/100;
	int64_t NROW = 1;
	int64_t NITER = 10*1000*1000 / NROW;

	// boost::random::mt19937 rng((unsigned int64_t)time(0));
	boost::random::mt19937 rng((uint64_t)0);	
	boost::random::uniform_int_distribution<> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	for(int64_t i = 0; i < NFILL; ++i) h[randindex(rng)] = i;
	cout << "done fill " 
	     << (double)h.size()/1000000 << "M  "
	     << (double)h.bucket_count()/1000000000 * sizeof(typename Map::value_type) << "GB  " 
	     << (double)h.size() / h.bucket_count()
	     << endl;	     	     

	boost::timer::cpu_timer t;
	size_t count = 0;
	for(size_t i = 0; i < NITER; ++i){
		size_t ri = randindex(rng);
		for(size_t j = 0; j < NROW; ++j){
			if( h.find(ri+j) != h.end() ) count += h[ri+j][0];
		}
	}
	cout << "rate " << (double)t.elapsed().wall/NITER*NROW << "ns, nonsense: " << (double)count / NITER << endl;
}

template<class MAP1, class MAP2>
void test_map(
	MAP1 & m, std::string lm,
	MAP2 & n, std::string ln
){
	int64_t MAXIDX = 10*1000*1000;
	int64_t NFILL = MAXIDX/10;
	int64_t NROW = 1;
	int64_t NITER = 100*1000*1000 / NROW;

	// boost::random::mt19937 rng((unsigned int64_t)time(0));
	boost::random::mt19937 rng((uint64_t)0);	
	boost::random::uniform_int_distribution<> randindex(0,MAXIDX);

	for(int64_t i = 0; i < NFILL; ++i){
		size_t ri = randindex(rng);
		m[ri] = 1;
		n[ri] = 1;
	}

	std::vector<size_t> ridx;
	for(size_t i = 0; i < NITER; ++i)
		ridx.push_back( randindex(rng) );

	boost::timer::cpu_timer tm;
	size_t mcount = 0;
	BOOST_FOREACH(size_t ri,ridx){
		for(size_t j = 0; j < NROW; ++j){
			if( m.find(ri+j) != m.end() ) mcount += m[ri+j];
		}
	}
	cout << lm << " rate " << (double)tm.elapsed().wall/NITER*NROW << "ns" << endl;

	boost::timer::cpu_timer tn;
	size_t ncount = 0;
	BOOST_FOREACH(size_t ri,ridx){
		for(size_t j = 0; j < NROW; ++j){
			if( n.find(ri+j) != n.end() ) ncount += n[ri+j];
		}
	}
	cout << ln << " rate " << (double)tn.elapsed().wall/NITER*NROW << "ns" << endl;

	cout << mcount << endl;
	ASSERT_EQ( mcount, ncount );
}

// 100M 1/10 nrow 1
// dense rate 1.04204e+07 -6.16073
// std   rate 4.29674e+06 -18.0023
// boost rate 5.60315e+06 4.98979

// 100M 1/10 nrow 10
// dense rate 2.15978e+07 176.247
// boost rate 6.85248e+06 200.992

struct GoogleDense { template<class K, class V> struct apply { typedef google::dense_hash_map<K,V> type; }; };

TEST(test_hash, DISABLED_sparse_vs_dense){
	SegmentedMap<uint64_t,util::SimpleArray<8,double> ,GoogleDense,2> m;

	google::dense_hash_map<uint64_t,util::SimpleArray<8,double> > d;
	d.set_empty_key(std::numeric_limits<uint64_t>::max());

	google::sparse_hash_map<uint64_t,util::SimpleArray<8,double> > s;

	// test_map(n,"google_dense",m,"segment_gdh ");
	cout << "====================== DENSE ======================" << endl;
	test_map(d); d.clear();
	cout << "====================== SPARSE =====================" << endl;
	test_map(s); s.clear();
}




// TEST(test_hash,std_unordered_map){
// 	std::unordered_map<int64_t,int64_t> h;
// 	test_map(h);
// }

// TEST(test_hash,boost_unordered_map){
// 	boost::unordered_map<int64_t,int64_t> h;
// 	test_map(h);
// }

// TEST(test_hash,std_map){
// 	std::map<int64_t,int64_t> h;
// 	test_map(h);
// }


}}
