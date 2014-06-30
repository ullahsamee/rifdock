#ifndef INCLUDED_util_mpl_TypesHolders_HH
#define INCLUDED_util_mpl_TypesHolders_HH

#include <boost/mpl/empty.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <vector>
#include <functional>

namespace scheme {
namespace util {
namespace mpl {

///@brief utility for mapping a Type into an empty placeholder that can be optimized away
template<typename T> struct type2type { typedef T type; };

namespace types_holder_details {
    template<
        typename Types,
        bool Empty = false
    >
    struct TypesValuesBase :
        public TypesValuesBase<
                    typename boost::mpl::pop_front<Types>::type,
                    boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
               >
    {
        typedef TypesValuesBase<
                     typename boost::mpl::pop_front<Types>::type,
                     boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
                > Super;
        typedef typename boost::mpl::front<Types>::type Held_t;
    
        using Super::get_impl;
        Held_t val_;
    
        Held_t       & get_impl(type2type<Held_t>)       { return val_; }
        Held_t const & get_impl(type2type<Held_t>) const { return val_; }
    };
    template<typename Types>
    struct TypesValuesBase<Types,true> {
        void get_impl(){}
    };
}

///@brief holds instances of a list of Types and provides a generic accessor get<T>()
///@tparam Types a "list" of types to hold instances of, most boost::mpl type containers work
///@detail
template< typename Types >
struct TypesValues :
    public types_holder_details::TypesValuesBase< Types, boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value >
{
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Types> ));
    template<typename T> T       & get()       { return get_impl( type2type<T>() ); }
    template<typename T> T const & get() const { return get_impl( type2type<T>() ); }
};

namespace types_container_details {
    template<
        typename Types,
        template<typename T,typename Alloc> class Container = std::vector,
        template<typename T> class Alloc = std::allocator,
        bool Empty = false
    >
    struct TypesContainersBase :
        public TypesContainersBase<
                    typename boost::mpl::pop_front<Types>::type,
                    Container, Alloc,
                    boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
               >
    {
        typedef TypesContainersBase<
                     typename boost::mpl::pop_front<Types>::type,
                     Container, Alloc,
                     boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
                > Super;
        typedef typename boost::mpl::front<Types>::type Held_t;
        typedef Container<Held_t,Alloc<Held_t> > Container_t;
    
        using Super::get_impl;
        Container_t val_;
    
        Container_t       & get_impl(type2type<Held_t>)       { return val_; }
        Container_t const & get_impl(type2type<Held_t>) const { return val_; }
    };
    template<
        typename Types,
        template<typename T,typename Alloc> class Container,
        template<typename T> class Alloc
    >
    struct TypesContainersBase<Types,Container,Alloc,true> {
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
struct TypesContainers :
    public types_container_details::TypesContainersBase< Types, Container, Alloc, boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value >
{
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Types> ));
    template<typename T> Container<T,Alloc<T> >       & get()       { return get_impl( type2type<T>() ); }
    template<typename T> Container<T,Alloc<T> > const & get() const { return get_impl( type2type<T>() ); }
};


namespace types_container_details {
    template<
        typename Types,
        template<typename T, typename Compare, typename Alloc> class Container = std::vector,
        template<typename T> class Compare = std::less,
        template<typename T> class Alloc = std::allocator,        
        bool Empty = false
    >
    struct TypesContainersCompBase :
        public TypesContainersCompBase<
                    typename boost::mpl::pop_front<Types>::type,
                    Container, Compare, Alloc, 
                    boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
               >
    {
        typedef TypesContainersCompBase<
                     typename boost::mpl::pop_front<Types>::type,
                     Container, Compare, Alloc,
                     boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
                > Super;
        typedef typename boost::mpl::front<Types>::type Held_t;
        typedef Container<Held_t,Compare<Held_t>,Alloc<Held_t> > Container_t;
    
        using Super::get_impl;
        Container_t val_;
    
        Container_t       & get_impl(type2type<Held_t>)       { return val_; }
        Container_t const & get_impl(type2type<Held_t>) const { return val_; }
    };
    template<
        typename Types,
        template<typename T, typename Compare, typename Alloc> class Container,
        template<typename T> class Compare,
        template<typename T> class Alloc        
    >
    struct TypesContainersCompBase<Types,Container,Compare,Alloc,true> {
        void get_impl(){}
    };
}

///@brief holds Containers of a list of Types and provides a generic accessor get<T>()
///@tparam Types a "list" of types to hold instances of, most boost::mpl type containers work
///@detail
template<
    typename Types,
    template<typename T, typename Compare, typename Alloc> class Container = std::vector,
    template<typename T> class Compare = std::less,
    template<typename T> class Alloc = std::allocator
>
struct TypesContainersComp :
    public types_container_details::TypesContainersCompBase< Types, Container, Compare, Alloc, boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value >
{
    BOOST_MPL_ASSERT(( boost::mpl::is_sequence<Types> ));     
    template<typename T> Container<T,Compare<T>,Alloc<T> >       & get()       { return get_impl( type2type<T>() ); }
    template<typename T> Container<T,Compare<T>,Alloc<T> > const & get() const { return get_impl( type2type<T>() ); }
};



// namespace types_svec_details {
//     template<
//         typename Types,
//         bool Empty = false
//     >
//     struct TypesVecSBase :
//         public TypesVecSBase<
//                     typename boost::mpl::pop_front<Types>::type,
//                     boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
//                >
//     {
//         typedef TypesVecSBase<
//                      typename boost::mpl::pop_front<Types>::type,
//                      boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value
//                 > Super;
//         typedef typename boost::mpl::front<Types>::type Held_t;
//         typedef std::vector<Held_t> Container_t;
//         typedef std::pair<typename Container_t::iterator,typename Container_t::iterator> IterPair;
//         typedef std::pair<typename Container_t::const_iterator,typename Container_t::const_iterator> ConstIterPair;

//         using Super::get_impl;
//         using Super::add_impl;
//         using Super::size_impl;    

//         Container_t vals_;

//         IterPair      get_impl(type2type<Held_t>)       { return make_pair(vals_.begin(),vals_.end()); }
//         ConstIterPair get_impl(type2type<Held_t>) const { return make_pair(vals_.begin(),vals_.end()); }

//         void add_impl( Held_t const & val ) { vals_.push_back(val); }

//         size_t size_impl( type2type<Held_t> ) { return vals_.size(); }
//     };
//     template<typename Types>
//     struct TypesVecSBase<Types,true> {
//         void get_impl(){}
//         void add_impl(){}    
//         void size_impl(){}    
//     };
// }

// template< typename Types >
// struct TypesVecS :
//     public types_svec_details::TypesVecSBase< Types, boost::mpl::empty<typename boost::mpl::pop_front<Types>::type>::value >
// {
//     template<typename T> std::pair<typename std::vector<T>::iterator      ,typename std::vector<T>::iterator      > get()       { return get_impl( type2type<T>() ); }
//     template<typename T> std::pair<typename std::vector<T>::const_iterator,typename std::vector<T>::const_iterator> get() const { return get_impl( type2type<T>() ); }

//     template<typename T> void add(T const & val) { return add_impl( val ); }

//     template<typename T> size_t size() { return size_impl( type2type<T>() ); }
// };



}
}
}

#endif
