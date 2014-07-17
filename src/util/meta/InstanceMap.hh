#ifndef INCLUDED_util_meta_InstanceMap_HH
#define INCLUDED_util_meta_InstanceMap_HH

#include <boost/mpl/zip_view.hpp>
#include <boost/mpl/is_sequence.hpp>

#include <boost/fusion/include/as_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>
#include <boost/fusion/include/mpl.hpp>

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
    
    template<typename K>
    typename boost::fusion::result_of::value_at_key<type,K>::type & 
    get() { return boost::fusion::at_key<K>(*this); }
    
    template<typename K>
    typename boost::fusion::result_of::value_at_key<type,K>::type const & 
    get() const { return boost::fusion::at_key<K>(*this); }
    
};



}
}
}

#endif
