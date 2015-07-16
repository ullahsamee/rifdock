#ifndef INCLUDED_objective_storage_TwoBodyTable_HH
#define INCLUDED_objective_storage_TwoBodyTable_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/util/assert.hh"

#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

namespace scheme { namespace objective { namespace storage {

template< class _Data = float >
struct TwoBodyTable {
	typedef _Data Data;
	typedef boost::multi_array< Data, 2 > Array2D;
	size_t nres_, nrot_;
	Array2D onebody_;
	boost::multi_array< int, 2 > all2sel_, sel2all_;
	std::vector<int> nsel_;
	boost::multi_array< Array2D, 2 > twobody_;
	TwoBodyTable( size_t nres, size_t nrots ) :
		nres_(nres),
		nrot_(nrots),
		onebody_( boost::extents[nres][nrots] ),
		all2sel_( boost::extents[nres][nrots] ),
		sel2all_( boost::extents[nres][nrots] ), // some at the ends of dim2 will be -1	
		nsel_( nres, 0 ),
		twobody_ ( boost::extents[nres][nres] )
	{}
	Data const & onebody( int ires, int irot ) const { return onebody_[ires][irot]; }
	void bounds_check_1b( int ires, int irot ) const {
		ALWAYS_ASSERT_MSG( ires >= 0, "ires < 0!" );
		ALWAYS_ASSERT_MSG( irot >= 0, "irot < 0!" );		
		ALWAYS_ASSERT_MSG( ires < nres_, "iret >= nres!" );
		ALWAYS_ASSERT_MSG( irot < nrot_, "irot >= nrot!" );				
	}
	Data const & onebody_at( int ires, int irot ) const {
		bounds_check_1b(ires,irot);
		return onebody_[ires][irot];
	}
	void set_onebody( int ires, int irot, Data const & val ){
		bounds_check_1b(ires,irot);
		onebody_[ires][irot] = val; 
	}
	Data twobody( int ires, int jres, int irot, int jrot ) const {
		if( twobody_[ires][jres].num_elements() ){
			int irotsel = all2sel_[ires][irot];
			int jrotsel = all2sel_[jres][jrot];
			if( irotsel < 0 || jrotsel < 0 ){
				return Data(9e9); // at least one of the onebodies is bad
			}
			return twobody_[ ires ][ jres ][ irotsel ][ jrotsel ];
		} else {
			return Data(0.0);
		}
	}
	// assumes onebody energies have been filled in at this point!
	void init_onebody_filter( float thresh ){
		for( int ires = 0; ires < nres_; ++ires ){
			int i = 0;
			for( int irot = 0; irot < nrot_; ++irot ){
				if( onebody_[ires][irot] <= thresh ){
					all2sel_[ires][irot] = i;
					sel2all_[ires][i] = irot;
					++i;
				} else {
					all2sel_[ires][irot] = -1;
				}
				nsel_[ires] = i;
				for( int j = i; j < nrot_; ++j ){
					sel2all_[ires][j] = -1; // past end of selected
				}
			}
		}
	}
	void init_twobody( int ires, int jres ){
		twobody_[ires][jres].resize( boost::extents[ nsel_[ires] ][ nsel_[jres] ] );
	}
	void clear_twobody( int ires, int jres ){
		twobody_[ires][jres].resize( boost::extents[0][0] );
	}

};

}}}

#endif
