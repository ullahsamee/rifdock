#include <gtest/gtest.h>
#include <util/mpl/TypesHolders.hh>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/begin.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/iterator_range.hpp>
#include <boost/mpl/for_each.hpp>
#include <vector>
#include <list>
#include <boost/foreach.hpp>

namespace scheme {
namespace util {
namespace mpl {

using std::cout;
using std::endl;
using std::vector;

struct PrintType {
    template<typename T>
    void operator()(T const & t){
        cout << typeid(t).name() << endl;
    }
};

TEST(TypesValues,holds_multiple_types_mpl_vector){
    TypesValues<boost::mpl::vector<double,int,char> > holder;
    holder.get<double>() = 1.0;
    holder.get<int>() = 2;  
    holder.get<char>() = 'A';
    ASSERT_EQ( holder.get<double>() , 1.0 );
    ASSERT_EQ( holder.get<int>()    , 2   );
    ASSERT_EQ( holder.get<char>()   , 'A' );
    holder.get<double>() = 1.23456789;
    holder.get<int>()    = 7;
    holder.get<char>()   = 'B';
    ASSERT_EQ( holder.get<double>() , 1.23456789 );
    ASSERT_EQ( holder.get<int>()    , 7   );
    ASSERT_EQ( holder.get<char>()   , 'B' );

    // boost::mpl::for_each<boost::mpl::vector<double,int,char> >( PrintType() );
}

TEST(TypesValues,holds_multiple_types_mpl_set){
    typedef boost::mpl::set<double,int,char> Types;
    TypesValues<Types> holder;
    holder.get<double>() = 1.0;
    holder.get<int>() = 2;  
    holder.get<char>() = 'A';
    ASSERT_EQ( holder.get<double>() , 1.0 );
    ASSERT_EQ( holder.get<int>()    , 2   );
    ASSERT_EQ( holder.get<char>()   , 'A' );
    holder.get<double>() = 1.23456789;
    holder.get<int>()    = 7;
    holder.get<char>()   = 'B';
    ASSERT_EQ( holder.get<double>() , 1.23456789 );
    ASSERT_EQ( holder.get<int>()    , 7   );
    ASSERT_EQ( holder.get<char>()   , 'B' );
}
TEST(TypesValues,holds_multiple_types_mpl_veiw){
    typedef boost::mpl::set<double,int,char> Types;
    typedef boost::mpl::iterator_range<
        boost::mpl::begin<Types>::type,
        boost::mpl::end  <Types>::type > View;
    TypesValues<View> holder;
    holder.get<double>() = 1.0;
    holder.get<int>() = 2;  
    holder.get<char>() = 'A';
    ASSERT_EQ( holder.get<double>() , 1.0 );
    ASSERT_EQ( holder.get<int>()    , 2   );
    ASSERT_EQ( holder.get<char>()   , 'A' );
    holder.get<double>() = 1.23456789;
    holder.get<int>()    = 7;
    holder.get<char>()   = 'B';
    ASSERT_EQ( holder.get<double>() , 1.23456789 );
    ASSERT_EQ( holder.get<int>()    , 7   );
    ASSERT_EQ( holder.get<char>()   , 'B' );
}

TEST(TypesValues,holds_multiple_vector_types_mpl_list){
    TypesValues<
        boost::mpl::list<
            vector<double>,
            vector<int>,
            vector<char>
        > > holder;
    holder.get<vector<double> >().push_back(1.0);
    holder.get<vector<int   > >().push_back( 2 );  
    holder.get<vector<char  > >().push_back('A');
    ASSERT_EQ( holder.get<vector<double> >().size() , 1 );
    ASSERT_EQ( holder.get<vector<int   > >().size() , 1 );
    ASSERT_EQ( holder.get<vector<char  > >().size() , 1 );
    ASSERT_EQ( holder.get<vector<double> >()[0] , 1.0 );
    ASSERT_EQ( holder.get<vector<int   > >()[0] ,  2  );
    ASSERT_EQ( holder.get<vector<char  > >()[0] , 'A' );
    holder.get<vector<double> >().push_back(2.0);
    holder.get<vector<int   > >().push_back( 6 );  
    holder.get<vector<char  > >().push_back('B');
    ASSERT_EQ( holder.get<vector<double> >()[0] , 1.0 );
    ASSERT_EQ( holder.get<vector<int   > >()[0] ,  2  );
    ASSERT_EQ( holder.get<vector<char  > >()[0] , 'A' );
    ASSERT_EQ( holder.get<vector<double> >()[1] , 2.0 );
    ASSERT_EQ( holder.get<vector<int   > >()[1] ,  6  );
    ASSERT_EQ( holder.get<vector<char  > >()[1] , 'B' );
}

TEST(TypesContainers,holds_multiple_types_mpl_list){
    TypesContainers<
        boost::mpl::list<
            double,
            int,
            char
        >, std::vector > holder;
    holder.get<double>().push_back(1.0);
    holder.get<int   >().push_back( 2 );  
    holder.get<char  >().push_back('A');
    ASSERT_EQ( holder.get<double>().size() , 1 );
    ASSERT_EQ( holder.get<int   >().size() , 1 );
    ASSERT_EQ( holder.get<char  >().size() , 1 );
    ASSERT_EQ( holder.get<double>()[0] , 1.0 );
    ASSERT_EQ( holder.get<int   >()[0] ,  2  );
    ASSERT_EQ( holder.get<char  >()[0] , 'A' );
    holder.get<double>().push_back(2.0);
    holder.get<int   >().push_back( 6 );  
    holder.get<char  >().push_back('B');
    ASSERT_EQ( holder.get<double>().size() , 2 );
    ASSERT_EQ( holder.get<int   >().size() , 2 );
    ASSERT_EQ( holder.get<char  >().size() , 2 );
    ASSERT_EQ( holder.get<double>()[0] , 1.0 );
    ASSERT_EQ( holder.get<int   >()[0] ,  2  );
    ASSERT_EQ( holder.get<char  >()[0] , 'A' );
    ASSERT_EQ( holder.get<double>()[1] , 2.0 );
    ASSERT_EQ( holder.get<int   >()[1] ,  6  );
    ASSERT_EQ( holder.get<char  >()[1] , 'B' );
}

TEST(TypesContainers,holds_lists){
    TypesContainers<
        boost::mpl::vector<
            double,
            int,
            char
        >, std::list > holder;
    holder.get<double>().push_back(1.0);
    holder.get<int   >().push_back( 2 );  
    holder.get<char  >().push_back('A');
    ASSERT_EQ( holder.get<double>().size() , 1 );
    ASSERT_EQ( holder.get<int   >().size() , 1 );
    ASSERT_EQ( holder.get<char  >().size() , 1 );
    ASSERT_EQ( holder.get<double>().front() , 1.0 );
    ASSERT_EQ( holder.get<int   >().front() ,  2  );
    ASSERT_EQ( holder.get<char  >().front() , 'A' );
    holder.get<double>().push_back(2.0);
    holder.get<int   >().push_back( 6 );  
    holder.get<char  >().push_back('B');
    ASSERT_EQ( holder.get<double>().size() , 2 );
    ASSERT_EQ( holder.get<int   >().size() , 2 );
    ASSERT_EQ( holder.get<char  >().size() , 2 );
    ASSERT_EQ( *(holder.get<double>().begin()) , 1.0 );
    ASSERT_EQ( *(holder.get<int   >().begin()) ,  2  );
    ASSERT_EQ( *(holder.get<char  >().begin()) , 'A' );
    ASSERT_EQ( *(++holder.get<double>().begin()) , 2.0 );
    ASSERT_EQ( *(++holder.get<int   >().begin()) ,  6  );
    ASSERT_EQ( *(++holder.get<char  >().begin()) , 'B' );
}


}
}
}
