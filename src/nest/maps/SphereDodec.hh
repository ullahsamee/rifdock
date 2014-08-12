#ifndef INCLUDED_scheme_nest_maps_SphereDodec_HH
#define INCLUDED_scheme_nest_maps_SphereDodec_HH

#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace maps {

	using namespace Eigen;

	namespace dodec_data {
		double * get_dodec_cell_centers(){
			static double p = (1.0+sqrt(5.0)) / 2.0;
			static double dodec_cell_centers[12*3] = {
				+1,+p, 0,
				+1,-p, 0,
				-1,+p, 0,
				-1,-p, 0,
				+p, 0,+1,
				-p, 0,+1,
				+p, 0,-1,
				-p, 0,-1,
				 0,+1,+p,
				 0,+1,-p,
				 0,-1,+p,
				 0,-1,-p,
			};
			return dodec_cell_centers;
		}
	}


	template<
		int DIM,
		class Value=Eigen::Matrix<double,3,1>,
		class Index=size_t,
		class Float=double
	>
	struct SphereDodec {
		BOOST_STATIC_ASSERT_MSG(DIM==2,"SphereDodec DIM must be == 2");

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
			Index resl,
			Value & value
		) const {
			Map< Matrix<double,3,12> > const cell_centers( dodec_data::get_dodec_cell_centers() );
			std::cout << cell_centers << std::endl;
			return false;
		}

		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Index resl,			
			Params & params,
			Index & cell_index
		) const ;

		///@brief get parameter space repr of Value for particular cell
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Index resl,
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
		Index num_cells() const ;
	};



}
}
}

#endif
