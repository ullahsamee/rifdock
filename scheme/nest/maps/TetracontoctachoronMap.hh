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

		Index nside_;
		Float one_over_nside_;

		TetracontoctachoronMap() : nside_(1), one_over_nside_(1.0) {}
		TetracontoctachoronMap(Index nside) : nside_(nside), one_over_nside_(1.0/nside) {}

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

			Params p = params * one_over_nside_;

			Index h48_cell_index = cell_index / (nside_*nside_*nside_);
			cell_index = cell_index % (nside_*nside_*nside_);
			p[0] += one_over_nside_ * (Float)( cell_index                   % nside_ );
			p[1] += one_over_nside_ * (Float)( cell_index /  nside_         % nside_ );
			p[2] += one_over_nside_ * (Float)( cell_index / (nside_*nside_) % nside_ );

			assert( p[0] >= 0.0 && p[0] <= 1.0 );
			assert( p[1] >= 0.0 && p[1] <= 1.0 );
			assert( p[2] >= 0.0 && p[2] <= 1.0 );						

			// std::cout << cell_index << " " << p << " " << p << std::endl;			
			// static int count = 0; if( ++count > 30 ) std::exit(-1);

			p = w*(p-0.5);

			// if( resl > 3 ){
				Float corner_dist = fabs(p[0])+fabs(p[1])+fabs(p[2]);
				Float delta = sqrt(w*w*w) / (Float)(1<<resl);
			// // static int count = 0;
   			// //          std::cout << corner_dist << "    " << p << " " << p << std::endl;
		   //          if(++count > 100) std::exit(-1);
			 	if( corner_dist - delta > 0.925 ) return false;
			// }

			Eigen::Quaternion<Float> q( sqrt(1.0-p.squaredNorm()), p[0], p[1], p[2] );
			assert( fabs(q.squaredNorm()-1.0) < 0.000001 );

			q = hbt24_cellcen<Float>( h48_cell_index ) * q;

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

			Index h48_cell_index;
			numeric::get_cell_48cell_half( q.coeffs(), h48_cell_index );

			q = hbt24_cellcen<Float>( h48_cell_index ).inverse() * q;

			q = to_half_cell(q);

			// // cout << q.w() << endl;
			// assert( q.w() > 0.7 );

			// // cout << "    get q0 " << q.coeffs().transpose() << endl;

			// // cout << "      get p  " << q.x() << " " << q.y() << " " << q.z() << endl;

			params[0] = q.x()/cell_width<Float>() + 0.5;
			params[1] = q.y()/cell_width<Float>() + 0.5;			
			params[2] = q.z()/cell_width<Float>() + 0.5;

			Indices ci = params * nside_;
			cell_index = ci[0] + ci[1]*nside_ + ci[2]*nside_*nside_;

			params = params * nside_ - ci.template cast<Float>();
			assert( params[0] >= 0.0 && params[0] <= 1.0 );
			assert( params[1] >= 0.0 && params[1] <= 1.0 );
			assert( params[2] >= 0.0 && params[2] <= 1.0 );						


			cell_index += h48_cell_index * nside_*nside_*nside_;

			// cout << "        get p0 " << params << endl;

			return true;
		}

		///@brief get parameter space repr of Value for particular cell
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Params & params
		) const {
			std::cerr << "Not Implemented" << std::endl;
			std::exit(-1);			
		}

		///@brief return the cell_index of neighboring cells within radius of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float radius, OutIter out) const {
			std::cerr << "Not Implemented" << std::endl;
			std::exit(-1);
		}

		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const {
			BOOST_VERIFY( resl < 6 );
			if( resl == 0 ){
				static Float const covrad[25] = {
					 62.76235 ,
					 38.63604 ,
					 26.71264 ,
					 20.62393 ,
					 17.02567 ,
					 14.25487 ,
					 12.42992 ,
					 11.02897 ,
					  9.62588 ,
					  8.70544 ,
					  7.82964 ,
					  7.28521 ,
					  6.62071 ,
					  6.13243 ,
					  5.81918 ,
					  5.44871 ,
					  5.14951 ,
					  4.82331 ,
					  4.52938 ,
					  4.31905 ,
					  4.07469 ,
					  3.93772 ,
					  3.77275 ,
					  3.64786 ,
					  3.44081 
				};
				return covrad[nside_-1];
			}
			assert(nside_==1);
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
		Index num_cells() const { return 24*nside_*nside_*nside_; }

	};



}}}

#endif
