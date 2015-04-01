#ifndef INCLUDED_scheme_nest_maps_OriTransMap_HH
#define INCLUDED_scheme_nest_maps_OriTransMap_HH

#include "scheme/nest/pmap/ScaleMap.hh"
#include "scheme/nest/pmap/TetracontoctachoronMap.hh"

#include "scheme/util/SimpleArray.hh"

#include <Eigen/Dense>

#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme { namespace nest { namespace pmap {

	template<
		int DIM=6,
		class Value=Eigen::Transform<double,3,Eigen::AffineCompact>,
		class Index=uint64_t,
		class Float=double
	>
	struct OriTransMap {
		BOOST_STATIC_ASSERT_MSG(DIM==6,"RotTransMap DIM must be == 6");

		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef util::SimpleArray<DIM,Index> Indices;
		typedef util::SimpleArray<DIM,Float> Params;
		typedef util::SimpleArray<3,Float> P3;
		typedef Eigen::Matrix<Float,3,1> V;
		typedef Eigen::Matrix<Float,3,3> M;


		typedef scheme::nest::pmap::TetracontoctachoronMap<3,M,Index,Float> OriMap;
		typedef scheme::nest::pmap::              ScaleMap<3,V,Index,Float> TransMap;

		OriMap ori_map_;
		TransMap trans_map_;

		template< class P, class I >
		OriTransMap(
			Float rot_resl_deg,
			P const & lb,
			P const & ub,
			I const & bs
		) : 
			ori_map_( OriMap::get_nside_for_rot_resl_deg( rot_resl_deg ) ),
			trans_map_( lb, ub, bs )
		{
			cout << "OriMap: TetracontoctachoronMap, " 
			     << rot_resl_deg << " nside: " << ori_map_.nside_ << ", covrad: " 
		    	 << ori_map_.bin_circumradius(0)*180.0/M_PI << ", size: " << ori_map_.num_cells() << endl;
		}


		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Index resl,
			Value & value
		) const {
			Index const ncori  = ori_map_.num_cells();
			Index const cori   = cell_index % ncori;
			Index const ctrans = cell_index / ncori;
			P3 pori  ( params, 0 ); // offset 0
			P3 ptrans( params, 3 ); // offset 3
			M m;
			V v;
			bool valid = ori_map_.params_to_value( pori  , cori  , resl, m );
			valid   &= trans_map_.params_to_value( ptrans, ctrans, resl, v );
			if( !valid ) return false;
			value = Value( m );
			value.translation()[0] = v[0];
			value.translation()[1] = v[1];
			value.translation()[2] = v[2];
			return true;
		}


		///@brief cell size
		Index num_cells() const { return ori_map_.num_cells() * trans_map_.num_cells(); }

	};



}}}

#endif
