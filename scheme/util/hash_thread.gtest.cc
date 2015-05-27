#ifdef CXX11

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

#include <thread>
#include <mutex>


namespace scheme { namespace util { namespace test_hash_thread {

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
void fill_map(Map & h, int64_t MAXIDX, int64_t sparsity=100ll ){
	int64_t NFILL = MAXIDX/sparsity;

	boost::random::mt19937 rng((uint64_t)0);	
	boost::random::uniform_int_distribution<int64_t> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	boost::timer::cpu_timer t;
	for(int64_t i = 0; i < NFILL; ++i) h[randindex(rng)] = i;
	cout << "done fill " 
		 << (double)h.size()/1000000 << "M  entries, "
		 << (double)h.bucket_count()/1000000000 * sizeof(typename Map::value_type) << "GB  total size, " 
		 << (double)h.size() / h.bucket_count() << " load, "
		 << (double)t.elapsed().wall/NFILL << "ns"
		 << endl;
}

template<class Map>
void test_map( Map * hp, double * runtime, int64_t MAXIDX, int64_t NITER ){
	Map & h(*hp);

	int64_t NROW=1;
	NITER /= NROW;

	static std::mutex m;

	// boost::random::mt19937 rng((unsigned int64_t)time(0));
	boost::random::mt19937 rng((uint64_t)0);	
	boost::random::uniform_int_distribution<int64_t> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	boost::timer::cpu_timer t;
	size_t count = 0;
	for(size_t i = 0; i < NITER; ++i){
		size_t ri = randindex(rng);
		for(size_t j = 0; j < NROW; ++j){
			typename Map::const_iterator iter = h.find(ri+j);
			count += iter==h.end() ? 0.0 : iter->second[0];
		}
	}
	double const time = (double)t.elapsed().wall;

	m.lock();
	*runtime += time;
	if( count == 0 ) cout << "ERROR!!" << endl;
	// cout << "rate " << (double)t.elapsed().wall/NITER/NROW << "ns, nonsense: " << (double)count << endl;
	m.unlock();
}

struct GoogleDense { template<class K, class V> struct apply { typedef google::dense_hash_map<K,V> type; }; };


TEST( test_hash_thread, google_hash_thread ){

 	int maxNthread = 12;
 	// SegmentedMap<uint64_t,util::SimpleArray<8,double> ,GoogleDense,2> m;

 	typedef google::dense_hash_map<uint64_t,util::SimpleArray<8,double> > D;
 	D d;
 	d.set_empty_key(std::numeric_limits<uint64_t>::max());

 	// typedef std::unordered_map<uint64_t,util::SimpleArray<8,double> > D;
 	// D d;

 	// typedef google::sparse_hash_map<uint64_t,util::SimpleArray<8,double> > D;
 	// D d;

 	// test_map(n,"google_dense",m,"segment_gdh ");
 	cout << "====================== HASH_TEST ======================" << endl;
 	int64_t MAXIDX = 4000ll*1000ll*1000ll;
 	int64_t NSAMP0 = 50ll*1000ll*1000ll;
 	fill_map( d, MAXIDX, 100 );

 	// double rt=0;
 	// test_map<D>( &d, &rt, MAXIDX, NSAMP0 ); 
 	// cout << "main thread: " << rt << "ns / lookup" << endl;

 	for(int Nthread = 1; Nthread <= maxNthread; ++Nthread){
 		int64_t NSAMP = NSAMP0 / Nthread;
 		// test_map(d); d.clear();
 		double runtime = 0.0;
 		std::vector<std::thread> t;
 		for(int i = 0; i < Nthread; ++i) t.push_back( std::thread( test_map<D>, &d, &runtime, MAXIDX, NSAMP ) );
 		for(int i = 0; i < Nthread; ++i) t[i].join();
 		runtime /= Nthread;
 		printf("nthread: %5i, %7.3fns / lookup, %7.3f runtime \n", Nthread, runtime / Nthread / NSAMP, runtime*1000000.0 );
 		std::cout.flush();
 	}	
 	d.clear();

 }


// 	// cout << "====================== SPARSE =====================" << endl;
// 	// test_map(s); s.clear();
// }

// ====================== dense_hash_map ======================
// done fill 39.8005M  9.66368GB  0.296537583.638ns
// main thread: 117.639
//                 dense       std    sparse
//                 9.7gb     3.8gb     3.0gb
// runtime     1 116.507   327.728   275.698
// runtime     2  47.328   158.355   115.406
// runtime     3  33.946   109.843    77.206
// runtime     4  24.162    86.198    58.460
// runtime     5  18.175    68.707    45.108
// runtime     6  16.432    57.417    39.264
// runtime     7  16.162    49.047    34.603
// runtime     8  23.386    46.593    34.271
// runtime     9  20.376    39.326    29.788
// runtime    10  15.526    32.694    27.142
// runtime    11  14.557    29.697    26.838
// runtime    12  14.538    29.900    26.953
// runtime    13  15.157    30.562    26.904
// runtime    14  14.894    32.081    25.822
// runtime    15  15.061    30.378    26.667
// runtime    16  15.734    31.775    27.109
// runtime    17  14.636    30.269    26.191
// runtime    18  14.716    29.658    25.827
// runtime    19  15.201    31.496    25.571
// runtime    20  14.980    34.030    25.837
// runtime    21  14.604    29.368    25.897
// runtime    22  14.157    28.862    25.810
// runtime    23  15.036    30.119    25.423
// runtime    24  14.653    31.193    26.734
// runtime    25  14.665    32.219    29.968
// runtime    26  14.626    31.588    30.335
// runtime    27  14.553    32.113    25.477
// runtime    28  14.699    34.224    26.479
// runtime    29  14.010    29.427    30.412
// runtime    30  14.786    29.614    27.794
// runtime    31  14.423    30.900    28.131
// runtime    32  14.913    32.057    27.559

// ====================== std::unordered_map ======================
// done fill 39.8005M  3.79296GB  0.755515664.705ns
// main thread: 315.627
// runtime     1 327.728
// runtime     2 158.355
// runtime     3 109.843
// runtime     4  86.198
// runtime     5  68.707
// runtime     6  57.417
// runtime     7  49.047
// runtime     8  46.593
// runtime     9  39.326
// runtime    10  32.694
// runtime    11  29.697
// runtime    12  29.900
// runtime    13  30.562
// runtime    14  32.081
// runtime    15  30.378
// runtime    16  31.775
// runtime    17  30.269
// runtime    18  29.658
// runtime    19  31.496
// runtime    20  34.030
// runtime    21  29.368
// runtime    22  28.862
// runtime    23  30.119
// runtime    24  31.193
// runtime    25  32.219
// runtime    26  31.588
// runtime    27  32.113
// runtime    28  34.224
// runtime    29  29.427
// runtime    30  29.614
// runtime    31  30.900
// runtime    32  32.057


// ====== sparsehash ACTUAL MEM ~3.0gb
// done fill 39.8005M  4.83184GB  0.5930741179.81ns
// main thread: 279.367
// runtime     1 275.698
// runtime     2 115.406
// runtime     3  77.206
// runtime     4  58.460
// runtime     5  45.108
// runtime     6  39.264
// runtime     7  34.603
// runtime     8  34.271
// runtime     9  29.788
// runtime    10  27.142
// runtime    11  26.838
// runtime    12  26.953
// runtime    13  26.904
// runtime    14  25.822
// runtime    15  26.667
// runtime    16  27.109
// runtime    17  26.191
// runtime    18  25.827
// runtime    19  25.571
// runtime    20  25.837
// runtime    21  25.897
// runtime    22  25.810
// runtime    23  25.423
// runtime    24  26.734
// runtime    25  29.968
// runtime    26  30.335
// runtime    27  25.477
// runtime    28  26.479
// runtime    29  30.412
// runtime    30  27.794
// runtime    31  28.131
// runtime    32  27.559

// dense_hash_map
// runtime     1 186.508 298.519
// runtime     2  85.484 141.353
// runtime     3  56.238  97.312
// runtime     4  43.064  74.743
// runtime     5  36.117  59.752
// runtime     6  33.066  53.263
// runtime     7  30.003  47.826
// runtime     8  32.299  45.008

// std::unordered_map
// runtime     1 298.519
// runtime     2 141.353
// runtime     3  97.312
// runtime     4  74.743
// runtime     5  59.752
// runtime     6  53.263
// runtime     7  47.826
// runtime     8  45.008
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


}}}

#endif
