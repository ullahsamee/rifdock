#ifndef INCLUDED_util_mpl_TypesHolders_HH
#define INCLUDED_util_mpl_TypesHolders_HH

#include <boost/mpl/deref.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/next.hpp>
#include <vector>

namespace scheme {
namespace util {
namespace mpl {

///@brief utility for mapping a Type into an empty placeholder that can be optimized away
template<typename T> struct type2type { typedef T type; };

namespace types_holder_details {
    using namespace boost::mpl;
    template< typename Beg, typename End >
    struct TypesValuesBase :
        public TypesValuesBase<typename next<Beg>::type, End >
    {
        typedef TypesValuesBase< typename next<Beg>::type, End > Super;
        typedef typename boost::mpl::deref<Beg>::type Held_t;
    
        using Super::get_impl;
        Held_t val_;
    
        Held_t       & get_impl(type2type<Held_t>)       { return val_; }
        Held_t const & get_impl(type2type<Held_t>) const { return val_; }
    };
    template<typename End>
    struct TypesValuesBase<End,End> {
        void get_impl(){}
    };
}

///@brief holds instances of a list of Types and provides a generic accessor get<T>()
///@tparam Types a "list" of types to hold instances of, most boost::mpl type containers work
///@detail
template< typename Types >
struct TypesValues : public types_holder_details::TypesValuesBase< 
                                typename boost::mpl::begin<Types>::type ,
                                typename boost::mpl::end  <Types>::type >
{
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Types> ));
    template<typename T> T       & get()       { return this->get_impl( type2type<T>() ); }
    template<typename T> T const & get() const { return this->get_impl( type2type<T>() ); }
};

namespace types_container_details {
    using namespace boost::mpl;
    template<
        typename Beg, typename End,
        template<typename T,typename Alloc> class Container = std::vector,
        template<typename T> class Alloc = std::allocator
    >
    struct TypesContainersBase : public TypesContainersBase< typename next<Beg>::type, End, Container, Alloc >
    {
        typedef TypesContainersBase< typename next<Beg>::type, End , Container, Alloc > Super;
        typedef typename deref<Beg>::type Held_t;
        typedef Container<Held_t,Alloc<Held_t> > Container_t;
    
        using Super::get_impl;
        Container_t val_;
    
        Container_t       & get_impl(type2type<Held_t>)       { return val_; }
        Container_t const & get_impl(type2type<Held_t>) const { return val_; }
    };
    template<
        typename End,
        template<typename T,typename Alloc> class Container,
        template<typename T> class Alloc
    >
    struct TypesContainersBase<End,End,Container,Alloc> {
        void get_impl(){}
    };
}

///@brief holds Containers of a list of Types and provides a generic accessor get<T>()
///@tparam Types a "list" of types to hold instances of, most boost::mpl type containers work
///@detail
template<
    typename Types,
    template<typename T,typename Alloc> class Container = std::vector,
    template<typename T> class Alloc = std::allocator
>
struct TypesContainers : public types_container_details::TypesContainersBase<
                                        typename boost::mpl::begin<Types>::type ,
                                        typename boost::mpl::end  <Types>::type ,
                                        Container, Alloc >
{
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Types> ));
    template<typename T> Container<T,Alloc<T> >       & get()       { return this->get_impl( type2type<T>() ); }
    template<typename T> Container<T,Alloc<T> > const & get() const { return this->get_impl( type2type<T>() ); }
};




}
}
}

#endif
