#ifndef INCLUDED_util_meta_InstanceMap_HH
#define INCLUDED_util_meta_InstanceMap_HH

#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <functional>
#include <boost/fusion/include/as_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/for_each.hpp>

namespace scheme {
namespace util {
namespace meta {

///@brief meta-container holding instances for any sequence of types
///@tparam Keys sequence of Key types
///@tparam Arg2 sequence of Value types OR metafunction class OR placeholder expression
///@detail if Arg2 is a metafunction class, the values that func applied to the Keys
template< typename Keys, typename Arg2 = Keys >
struct InstanceMap :
 boost::fusion::result_of::as_map<
    typename boost::mpl::transform<
            Keys
        ,   typename boost::mpl::eval_if< 
                    boost::mpl::is_sequence<Arg2>
                ,   Arg2
                ,   boost::mpl::transform<Keys,Arg2>
                >::type
        ,   boost::fusion::pair<boost::mpl::_1,boost::mpl::_2>
        >::type
    >::type
{
    typedef Keys KeyTypes;
    typedef typename boost::fusion::result_of::as_map<
        typename boost::mpl::transform<
            Keys
        ,   typename boost::mpl::eval_if< 
                    boost::mpl::is_sequence<Arg2>
                ,   Arg2
                ,   boost::mpl::transform<Keys,Arg2>
            >::type
        ,   boost::fusion::pair<boost::mpl::_1,boost::mpl::_2>
        >::type
    >::type type;
    
    InstanceMap(){}

    template<class A>
    InstanceMap(A const&a) : type(a) {}
    template<class A,class B>
    InstanceMap(A const&a,B const&b) : type(a,b) {}
    template<class A,class B,class C>
    InstanceMap(A const&a,B const&b,C const&c) : type(a,b,c) {}
    template<class A,class B,class C,class D>
    InstanceMap(A const&a,B const&b,C const&c,D const&d) : type(a,b,c,d) {}
    template<class A,class B,class C,class D,class E>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e) : type(a,b,c,d,e) {}
    template<class A,class B,class C,class D,class E,class F>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e,F const&f) : type(a,b,c,d,e,f) {}
    template<class A,class B,class C,class D,class E,class F,class G>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e,F const&f,G const&g) : type(a,b,c,d,e,f,g) {}
    template<class A,class B,class C,class D,class E,class F,class G,class H>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e,F const&f,G const&g,H const&h) : type(a,b,c,d,e,f,g,h) {}
    template<class A,class B,class C,class D,class E,class F,class G,class H,class I>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e,F const&f,G const&g,H const&h,I const&i) : type(a,b,c,d,e,f,g,h,i) {}
    template<class A,class B,class C,class D,class E,class F,class G,class H,class I,class J>
    InstanceMap(A const&a,B const&b,C const&c,D const&d,E const&e,F const&f,G const&g,H const&h,I const&i,J const&j) : type(a,b,c,d,e,f,g,h,i,j) {}

    template<typename K>
    typename boost::fusion::result_of::value_at_key<type,K>::type & 
    get() { return boost::fusion::at_key<K>(*this); }
    
    template<typename K>
    typename boost::fusion::result_of::value_at_key<type,K>::type const & 
    get() const { return boost::fusion::at_key<K>(*this); }
    
};

namespace impl {
    template<class Float> struct SUM { 
        Float & sum; SUM(Float & s):sum(s){}
        template<class T> void operator()(T const & x) const { sum += x.second; }
    };
    template<class Float> struct SETVAL { Float val; template<class T> void operator()(T & x) const { x.second = val; } };
    template<class T> struct ADD { T & sink; ADD(T & s) : sink(s) {}
        template<class X> void operator()(X const & x) const {
            sink.template get<typename X::first_type>() += x.second;
        }
    };
    template<class T> struct MUL { T & sink; MUL(T & s) : sink(s) {}
        template<class X> void operator()(X const & x) const {
            sink.template get<typename X::first_type>() *= x.second;
        }
    };
    template<class T,class OP> struct BINARY_OP_EQUALS { T & sink; BINARY_OP_EQUALS(T & s) : sink(s) {}
        template<class X> void operator()(X const & x) const {
            sink.template get<typename X::first_type>() = OP()(sink.template get<typename X::first_type>(),x.second);
        }
    };
}


template<class T> struct is_InstanceMap : boost::mpl::false_ {};
template<class A,class B> struct is_InstanceMap< InstanceMap<A,B> > : boost::mpl::true_ {};


/// values must be convertable into Float
template< typename Keys, typename Arg2 = Keys, class Float = double >
struct NumericInstanceMap : InstanceMap<Keys,Arg2> {
    typedef NumericInstanceMap<Keys,Arg2,Float> THIS;
    typedef InstanceMap<Keys,Arg2> BASE;
    typedef Float F;
    NumericInstanceMap(F f=0){ setall(f); }
    NumericInstanceMap(F a,F b) : BASE(a,b) {}
    NumericInstanceMap(F a,F b,F c) : BASE(a,b,c) {}
    NumericInstanceMap(F a,F b,F c,F d) : BASE(a,b,c,d) {}
    NumericInstanceMap(F a,F b,F c,F d,F e) : BASE(a,b,c,d,e) {}
    NumericInstanceMap(F a,F b,F c,F d,F e,F f) : BASE(a,b,c,d,e,f) {}
    NumericInstanceMap(F a,F b,F c,F d,F e,F f,F g) : BASE(a,b,c,d,e,f,g) {}
    NumericInstanceMap(F a,F b,F c,F d,F e,F f,F g,F h) : BASE(a,b,c,d,e,f,g,h) {}
    NumericInstanceMap(F a,F b,F c,F d,F e,F f,F g,F h,F i) : BASE(a,b,c,d,e,f,g,h,i) {}
    NumericInstanceMap(F a,F b,F c,F d,F e,F f,F g,F h,F i,F j) : BASE(a,b,c,d,e,f,g,h,i,j) {}
    void setall(Float val){
        impl::SETVAL<Float> set;
        set.val = val;
        boost::fusion::for_each( *this, set );
    }
    Float sum() const {
        Float sum=0;
        impl::SUM<Float> s(sum);
        boost::fusion::for_each( *this, s );
        return sum;
    }
    operator Float() const { return sum(); }
    void operator+=(THIS const & o){
        impl::BINARY_OP_EQUALS< THIS, std::plus<Float> > add(*this);
        boost::fusion::for_each(o,add);
    }
    void operator-=(THIS const & o){
        impl::BINARY_OP_EQUALS< THIS, std::minus<Float> > add(*this);
        boost::fusion::for_each(o,add);
    }
    void operator/=(THIS const & o){
        impl::BINARY_OP_EQUALS< THIS, std::divides<Float> > add(*this);
        boost::fusion::for_each(o,add);
    }
    void operator*=(THIS const & o){
        impl::BINARY_OP_EQUALS< THIS, std::multiplies<Float> > add(*this);
        boost::fusion::for_each(o,add);
    }
};

template<class A, class B, class C>
NumericInstanceMap<A,B,C> operator*(NumericInstanceMap<A,B,C> const & a, NumericInstanceMap<A,B,C> const & b){
    NumericInstanceMap<A,B,C> result = a;
    result *= b;
    return result;
}
template<class A, class B, class C>
NumericInstanceMap<A,B,C> operator+(NumericInstanceMap<A,B,C> const & a, NumericInstanceMap<A,B,C> const & b){
    NumericInstanceMap<A,B,C> result = a;
    result += b;
    return result;
}
template<class A, class B, class C>
NumericInstanceMap<A,B,C> operator-(NumericInstanceMap<A,B,C> const & a, NumericInstanceMap<A,B,C> const & b){
    NumericInstanceMap<A,B,C> result = a;
    result -= b;
    return result;
}
template<class A, class B, class C>
NumericInstanceMap<A,B,C> operator/(NumericInstanceMap<A,B,C> const & a, NumericInstanceMap<A,B,C> const & b){
    NumericInstanceMap<A,B,C> result = a;
    result /= b;
    return result;
}


}
}
}

#endif
