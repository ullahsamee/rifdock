#ifndef INCLUDED_scheme_util_SimpleArray_HH
#define INCLUDED_scheme_util_SimpleArray_HH

#include <types.hh>
#include <cmath>
#include <boost/assert.hpp>
#include <boost/serialization/access.hpp>

namespace scheme {
namespace util {

//TODO: finish this and replace Eigen dependency in NEST
namespace impl {
	struct NoInit {};
}

///@brief minimal fixed size array with element-wise operations
///@note used this instead of Eigen Array in NEST to speed compilation by 20%-50%
template<int N, class F=double>
struct SimpleArray {
	typedef F value_type;
	typedef F* iterator;
	typedef F const* const_iterator;
	typedef F & reference;
	typedef F const & const_reference;
	typedef size_t difference_type;
	typedef size_t size_type;
	F D[N];
	// SimpleArray() { for(size_t i = 0; i < N; ++i) D[i]=0; }
	SimpleArray(){}
	// explicit SimpleArray(F const* fp){ for(size_t i = 0; i < N; ++i) D[i] = fp[i]; }
	SimpleArray(F a){ fill(a); }
	SimpleArray(F a, F b                                        ){ D[0]=a;D[1]=b; }
	SimpleArray(F a, F b, F c                                   ){ D[0]=a;D[1]=b;D[2]=c; }
	SimpleArray(F a, F b, F c, F d                              ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d; }
	SimpleArray(F a, F b, F c, F d, F e                         ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e; }
	SimpleArray(F a, F b, F c, F d, F e, F f                    ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e;D[5]=f; }
	SimpleArray(F a, F b, F c, F d, F e, F f, F g               ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e;D[5]=f;D[6]=g; }
	SimpleArray(F a, F b, F c, F d, F e, F f, F g, F h          ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e;D[5]=f;D[6]=g;D[7]=h; }
	SimpleArray(F a, F b, F c, F d, F e, F f, F g, F h, F i     ){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e;D[5]=f;D[6]=g;D[7]=h;D[8]=i; }
	SimpleArray(F a, F b, F c, F d, F e, F f, F g, F h, F i, F j){ D[0]=a;D[1]=b;D[2]=c;D[3]=d;D[4]=e;D[5]=f;D[6]=g;D[7]=h;D[8]=i;D[9]=j; }
	F       & operator[](size_t i)       { return D[i]; }
	F const & operator[](size_t i) const { return D[i]; }
	F       & at(size_t i)       { BOOST_VERIFY(i < N); return D[i]; }
	F const & at(size_t i) const { BOOST_VERIFY(i < N); return D[i]; }
	template<class OF> SimpleArray<N,OF> cast() const {
		SimpleArray<N,OF> r; for(int i = 0; i < N; ++i) r[i] = (*this)[i]; return r;
	}
	SimpleArray<N,F> max(F b) const { SimpleArray<N,F> r(*this); for(int i = 0; i < N; ++i) r[i] = std::max(r[i],b); return r; }
	SimpleArray<N,F> min(F b) const { SimpleArray<N,F> r(*this); for(int i = 0; i < N; ++i) r[i] = std::min(r[i],b); return r; }
	SimpleArray<N,F> max(SimpleArray<N,F> const & b) const { SimpleArray<N,F> r(*this); for(int i = 0; i < N; ++i) r[i] = std::max(r[i],b[i]); return r; }
	SimpleArray<N,F> min(SimpleArray<N,F> const & b) const { SimpleArray<N,F> r(*this); for(int i = 0; i < N; ++i) r[i] = std::min(r[i],b[i]); return r; }
	F prod() const { F p=1; for(int i = 0; i < N; ++i) p *= D[i]; return p; }
	F sum () const { F p=0; for(int i = 0; i < N; ++i) p += D[i]; return p; }
	F prod(size_t l) const { F p=1; for(int i = 0; i < l; ++i) p *= D[i]; return p; }
	F sum (size_t l) const { F p=0; for(int i = 0; i < l; ++i) p += D[i]; return p; }
	bool operator==(SimpleArray<N,F> const & o) const {
		bool r = true; for(int i = 0; i < N; ++i) r &= D[i]==o.D[i]; return r;
	}
	F norm2() const { F n=0; for(int i = 0; i < N; ++i) n+=D[i]*D[i]; return n; }
	F norm() const { return std::sqrt(norm2()); }
	void fill(F v){ for(int i = 0; i < N; ++i) D[i]=v; }
	iterator begin() { return &D[0]; }
	iterator end  () { return &D[N]; }
	const_iterator begin() const { return &D[0]; }
	const_iterator end  () const { return &D[N]; }
	bool empty() const { return false; }
	size_type size() const { return N; }
	void swap(SimpleArray<N,F> & o){ for(int i = 0; i < N; ++i) std::swap(D[i],o.D[i]); }
	friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int ){
		for(size_t i = 0; i < N; ++i) ar & D[i];
	}
};
template<int N, class F>
std::ostream & operator<<(std::ostream & out,SimpleArray<N,F> const & a){
	for(int i = 0; i < N; ++i) out << a[i] << " ";
	return out;
}

template<int N, class F> SimpleArray<N,F> operator+( SimpleArray<N,F> const & a, SimpleArray<N,F> const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] += b[i]; return r; }
template<int N, class F> SimpleArray<N,F> operator-( SimpleArray<N,F> const & a, SimpleArray<N,F> const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] -= b[i]; return r; }
template<int N, class F> SimpleArray<N,F> operator*( SimpleArray<N,F> const & a, SimpleArray<N,F> const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] *= b[i]; return r; }
template<int N, class F> SimpleArray<N,F> operator/( SimpleArray<N,F> const & a, SimpleArray<N,F> const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] /= b[i];	return r; }
template<int N, class F> SimpleArray<N,F> operator+( SimpleArray<N,F> const & a, F const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] += b; return r; }
template<int N, class F> SimpleArray<N,F> operator-( SimpleArray<N,F> const & a, F const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] -= b; return r; }
template<int N, class F> SimpleArray<N,F> operator*( SimpleArray<N,F> const & a, F const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] *= b; return r; }
template<int N, class F> SimpleArray<N,F> operator/( SimpleArray<N,F> const & a, F const & b ){
	SimpleArray<N,F> r(a); for(int i = 0; i < N; ++i) r[i] /= b;	return r; }

}
}

#endif
