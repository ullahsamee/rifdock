#ifndef INCLUDED_scheme_nest_maps_SphereQuad_HH
#define INCLUDED_scheme_nest_maps_SphereQuad_HH

#include <nest/maps/parameter_maps.hh>
#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace maps {

	using namespace Eigen;
	using std::cout;
	using std::endl;

	namespace quadsphere_data {
		double * get_quadsphere_cells(){
			static double quadsphere_cells[9*6] = {
				 1, 0, 0,
	 			 0, 1, 0,
	 			 0, 0, 1,

				 1, 0, 0,
				 0, 0,-1,
				 0, 1, 0,

				 0, 0, 1,
				 0, 1, 0,
				-1, 0, 0,
				 
				 1, 0, 0,
				 0, 0, 1,
				 0,-1, 0,
				 
				 0, 0,-1,
				 0, 1, 0,
				 1, 0, 0,
				 
				 1, 0, 0,
				 0,-1, 0,
				 0, 0,-1,
			};
			return quadsphere_cells;
		}
	}


	template<
		int DIM,
		class Value=Eigen::Matrix<double,3,1>,
		class Index=size_t,
		class Float=double
	>
	struct SphereQuad {
		BOOST_STATIC_ASSERT_MSG(DIM==2,"SphereQuad DIM must be == 2");

		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<Float,DIM,1> Params;		

		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Value & value
		) const {
			Map< Matrix<double,3,3> > const rot_to_cell(
				const_cast<double*>( quadsphere_data::get_quadsphere_cells() + 9*cell_index ));
			Vector3d vec( params[0]*2-1 , params[1]*2-1 , 1 );
			vec = rot_to_cell * vec;
			vec = vec / vec.norm();
			value[0] = vec[0];
			value[1] = vec[1];
			value[2] = vec[2];						
			return true;
		}

		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			Matrix<double,3,1> tmpval( value[0], value[1], value[2] );
			tmpval = tmpval / tmpval.norm();
			// cout << tmpval.transpose() << endl;

			Map< Matrix<double,3,6>, 0, OuterStride<9> > const 
				centers( quadsphere_data::get_quadsphere_cells()+6 );
			// cout << centers << endl;
			Matrix<double,1,6> dotprods = centers.transpose() * tmpval;

			double highest = -9e9;
			for(size_t i = 0; i < 6; ++i){
				if( dotprods[i] > highest ){
					highest = dotprods[i];
					cell_index = i;
				}
			}
			Map< Matrix<double,3,3> > const rot_to_cell( 
				quadsphere_data::get_quadsphere_cells() + 9*cell_index );
			tmpval = rot_to_cell.transpose() * tmpval;
			params[0] = (tmpval[0]/highest+1.0)/2.0;
			params[1] = (tmpval[1]/highest+1.0)/2.0;
			// cout << "PRM " << params.transpose() << " " << highest << endl;
			assert( 0.0 <= params[0] && params[0] <= 1.0);
			assert( 0.0 <= params[1] && params[1] <= 1.0);			
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
			double const delta = 1.0/(double)(1ul<<resl);
			Vector3d pworst = Vector3d(1,1,1) - delta*Vector3d(2.0,0,0);
			Vector3d pNest0 = Vector3d(1,1,1) - delta*Vector3d(1.0,1.0,0);
			pworst = pworst / pworst.norm();
			pNest0 = pNest0 / pNest0.norm();
			return (pworst-pNest0).norm() * 1.001; // fudge factor
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
		Index cell_size() const { return 6; }
	};



}
}
}

#endif
