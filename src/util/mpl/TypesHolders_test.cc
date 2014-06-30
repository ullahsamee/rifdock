#include <gtest/gtest.h>
#include <util/mpl/TypesHolders.hh>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/list.hpp>
#include <vector>
#include <list>
#include <set>
#include <boost/foreach.hpp>

namespace scheme {
namespace util {
namespace mpl {

using std::cout;
using std::endl;
using std::vector;


TEST(TypesValues,holds_multiple_types_mplvector){
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

TEST(TypesContainers,holds_sets){
    TypesContainersComp<
        boost::mpl::vector<
            double,
            int,
            char
        >, std::set > holder;
    holder.get<double>().insert(1.0);
    holder.get<int   >().insert( 2 );  
    holder.get<char  >().insert('A');
    ASSERT_EQ( holder.get<double>().size() , 1 );
    ASSERT_EQ( holder.get<int   >().size() , 1 );
    ASSERT_EQ( holder.get<char  >().size() , 1 );
    ASSERT_EQ( *(holder.get<double>().begin()) , 1.0 );
    ASSERT_EQ( *(holder.get<int   >().begin()) ,  2  );
    ASSERT_EQ( *(holder.get<char  >().begin()) , 'A' );
    holder.get<double>().insert(2.0);
    holder.get<int   >().insert( 6 );  
    holder.get<char  >().insert('B');
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

// TEST(TypesContainer,holds_multiple_types){
//     TypesVecS<
//         boost::mpl::list<
//             double,
//             int,
//             char
//         > > holder;
//     ASSERT_EQ( holder.size<double>() , 0 );
//     ASSERT_EQ( holder.size<int   >() , 0 );
//     ASSERT_EQ( holder.size<char  >() , 0 );
//     holder.add<double>(1.0);
//     holder.add<int   >( 2 );  
//     holder.add<char  >('A');
//     ASSERT_EQ( holder.size<double>() , 1 );
//     ASSERT_EQ( holder.size<int   >() , 1 );
//     ASSERT_EQ( holder.size<char  >() , 1 );
//     std::pair<std::vector<double>::iterator,std::vector<double>::iterator> ip = holder.get<double>();
//     BOOST_FOREACH( double const & d, holder.get<double>() ){
//         cout << d << endl;
//     }
//     // ASSERT_EQ( holder.get<vector<double> >()[0] , 1.0 );
//     // ASSERT_EQ( holder.get<vector<int   > >()[0] ,  2  );
//     // ASSERT_EQ( holder.get<vector<char  > >()[0] , 'A' );
//     holder.add<double>(2.0);
//     holder.add<int   >( 6 );  
//     holder.add<char  >('B');
//     ASSERT_EQ( holder.size<double>() , 2 );
//     ASSERT_EQ( holder.size<int   >() , 2 );
//     ASSERT_EQ( holder.size<char  >() , 2 );
//     // ASSERT_EQ( holder.get<vector<double> >()[0] , 1.0 );
//     // ASSERT_EQ( holder.get<vector<int   > >()[0] ,  2  );
//     // ASSERT_EQ( holder.get<vector<char  > >()[0] , 'A' );
//     // ASSERT_EQ( holder.get<vector<double> >()[1] , 2.0 );
//     // ASSERT_EQ( holder.get<vector<int   > >()[1] ,  6  );
//     // ASSERT_EQ( holder.get<vector<char  > >()[1] , 'B' );

// }

}
}
}
