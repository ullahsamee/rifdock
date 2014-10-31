#ifndef INCLUDED_scheme_nest_maps_HecatonicosachoronMap_HH
#define INCLUDED_scheme_nest_maps_HecatonicosachoronMap_HH

#include "scheme/numeric/geom_4d.hh"

#include <Eigen/Dense>
#include "scheme/util/SimpleArray.hh"

#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme { namespace nest { namespace maps {

	using namespace Eigen;
	using std::cout;
	using std::endl;

	template<class Float> Float cell_width(){ return 0.766; }

	template<class Q> Q to_half_cell(Q const & q){ return q.w()>=0 ? q : Q(-q.w(),-q.x(),-q.y(),-q.z()); }

	template<class Float, class Index>
	Eigen::Map<Eigen::Quaternion<Float>const> hbt24_cellcen(Index const & i){
		// Float const * tmp = numeric::get_raw_48cell_half<Float>() + 4*i;
		// std::cout << "   raw hbt24_cellcen " << tmp[0] << " " << tmp[1] << " " << tmp[2] << " " << tmp[3] << std::endl;
		return Eigen::Map<Eigen::Quaternion<Float>const>(
		 numeric::get_raw_48cell_half<Float>() + 4*i );
	}

	template<
		int DIM=3,
		class Value=Eigen::Matrix3d,
		class Index=uint64_t,
		class Float=double
	>
	struct TetracontoctachoronMap {
		BOOST_STATIC_ASSERT_MSG(DIM==3,"TetracontoctachoronMap DIM must be == 3");

		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef util::SimpleArray<DIM,Index> Indices;
		typedef util::SimpleArray<DIM,Float> Params;		

		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Index resl,
			Value & value
		) const {
			// // cout << "        set p0 " << params << endl;
			Float const & w(cell_width<Float>());
			Params p = w*(params-0.5);


			if( resl > 3 ){
				Float corner_dist = fabs(p[0])+fabs(p[1])+fabs(p[2]);
				Float delta = sqrt(w*w*w) / (Float)(1<<resl);
			// // static int count = 0;
   			// //          std::cout << corner_dist << "    " << params << " " << p << std::endl;
		   //          if(++count > 100) std::exit(-1);
			 	if( corner_dist - delta > 0.925 ) return false;
			}

			Eigen::Quaternion<Float> q( sqrt(1.0-p.squaredNorm()), p[0], p[1], p[2] );
			assert( fabs(q.squaredNorm()-1.0) < 0.000001 );

			q = hbt24_cellcen<Float>(cell_index)*q;

			value = q.matrix();
			return true;
		}

		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Index /*resl*/,
			Params & params,
			Index & cell_index
		) const {
			Quaternion<Float> q(value);
			// q = to_half_cell(q);
			// // cout << "get q  " << q.coeffs().transpose() << endl;

			numeric::get_cell_48cell_half(q.coeffs(),cell_index);

			q = hbt24_cellcen<Float>(cell_index).inverse() * q;

			q = to_half_cell(q);

			// // cout << q.w() << endl;
			// assert( q.w() > 0.7 );

			// // cout << "    get q0 " << q.coeffs().transpose() << endl;

			// // cout << "      get p  " << q.x() << " " << q.y() << " " << q.z() << endl;

			params[0] = q.x()/cell_width<Float>() + 0.5;
			params[1] = q.y()/cell_width<Float>() + 0.5;			
			params[2] = q.z()/cell_width<Float>() + 0.5;

			// cout << "        get p0 " << params << endl;

			return true;
		}

		///@brief get parameter space repr of Value for particular cell
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Params & params
		) const ;

		///@brief return the cell_index of neighboring cells within radius of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float radius, OutIter out) const;

		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const {
			BOOST_VERIFY( resl < 6 );
			static Float const covrad[6] = {
				 62.76235, 
				 37.95720, 
				 20.53126, 
				 11.00814, 
				  5.31355, 
				  2.66953 
			};
			return covrad[resl];
		}

		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const {
			// double const delta = 1.0/(double)(1ul<<resl);
			// Vec pworst = Vec(1,1,1);// - delta*Vec(2.0,0,0);
			// Vec pNest0 = Vec(1,1,1) - delta*Vec(1.0,1.0,0);
			// cube_to_sphere(pworst);
			// cube_to_sphere(pNest0);
			// pworst.normalize();
			// pNest0.normalize();
			// return pworst.distance(pNest0) * 0.5999; // should be half of curcumradius based on geometry
		}

		///@brief cell size
		Index num_cells() const { return 24; }

	};



}}}

#endif
