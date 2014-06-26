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

	namespace sphere_quad_data {
		static double p = (1.0+sqrt(5.0)) / 2.0;
		static double const sphere_quad_cells[9*6] = {
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
			 0,-1,-0,
			 0, 0,-1,
		};
	}


	template<
		int DIM,
		class Value=Eigen::Matrix<double,3,1>,
		class Index=size_t,
		class Float=double
	>
	struct SphereQuad : public ParamMap<DIM,Value,Index,Float> {
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
				const_cast<double*>( sphere_quad_data::sphere_quad_cells + 9*cell_index ));
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
			cell_index = 1;
			Map< Matrix<double,3,3> > const rot_to_cell(
				const_cast<double*>( sphere_quad_data::sphere_quad_cells + 9*cell_index ));

			double test[18] = {
				1,2,3,
				4,5,6,
				7,8,9,
				10,11,12,
				13,14,15,
				16,17,18
			};
			Map< Matrix<double,3,2>, 0, OuterStride<9> > const centers(
				const_cast<double*>( test ));

			cout << "centers" << endl;
			cout << centers << endl;

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
		Float bin_circumradius(Index resl) const ;

		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const ;

		///@brief cell size
		Index cell_size() const ;
	};



}
}
}

#endif
