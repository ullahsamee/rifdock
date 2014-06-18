#ifndef INCLUDED_scheme_nest_parameter_maps_HH
#define INCLUDED_scheme_nest_parameter_maps_HH

#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <vector>

namespace scheme {
namespace nest {

	using namespace Eigen;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///	Parameter Mappings
	///////////////////////////////////////////////////////////////////////////////////////////////



	///@brief Parameter to Value Map Policy Class
	///@detail just copies the [0.0,1.0] hypercube coords to Value
	///        the first dimension will incremented by cell_index
	template< int DIM, class Value, class Index, class Float >
	struct UnitMap {
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<Float,DIM,1> Params;		
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Value & value
		) const {
			for(size_t i = 0; i < DIM; ++i) value[i] = params[i];
			value[0] += (Float)cell_index;
			return true;
		}
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			for(size_t i = 0; i < DIM; ++i) params[i] = value[i];
			cell_index = (Index)value[0];
			params[0] -= (Float)cell_index;
			return true;
		}
		Float covering_radius(Index resl) const { return 0.5/(Float)(1<<resl) * sqrt(DIM); }
		Float neighbor_radius(Index resl) const { return 1.5/(Float)(1<<resl); }
	 protected:
		~UnitMap(){}
	};
	


	/// @brief Parameter Mapping Policy for cartesian grids
	/// @note NEST cell_size MUST agree with cell_sizes_
	template< int DIM, class Value, class Index, class Float >
	struct ScaleMap {
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<Float,DIM,1> Params;		
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
	 private:
		///@brief lower bound on value space		
		Params lower_bound_;
		///@brief upper bound on value space base size 1
		Params upper_bound_;
		///@brief distributes cell_index accross dimensions
		Indices cell_sizes_;
		Indices cell_sizes_pref_sum_;
	 public:
		///@brief construct with default lb, ub, bs
		ScaleMap(){	cell_sizes_.fill(1); lower_bound_.fill(0); upper_bound_.fill(1); init(); }
		///@brief construct with default lb, ub
		ScaleMap(Indices const & bs) : cell_sizes_(bs) { lower_bound_.fill(0); upper_bound_.fill(1); init(); }
		///@brief construct with default bs
		ScaleMap(Params const & lb, Params const & ub) : lower_bound_(lb), upper_bound_(ub) { cell_sizes_.fill(1); init(); }
		///@brief construct with specified lb, ub and bs
		ScaleMap(Params const & lb, Params const & ub, Indices const & bs) : lower_bound_(lb), upper_bound_(ub), cell_sizes_(bs) { init(); }

		///@brief sets up cell_size_pref_sum
		void init(){
			for(size_t i = 0; i < DIM; ++i) 
				cell_sizes_pref_sum_[i] = cell_sizes_.head(i).prod();
		}

		///@brief sets value based on cell_index and parameters using geometric bounds
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Value & value
		) const {
			for(size_t i = 0; i < DIM; ++i){
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				Float bi = ( cell_index / cell_sizes_pref_sum_[i] ) % cell_sizes_[i];
				Float width = (upper_bound_[i]-lower_bound_[i])/(Float)cell_sizes_[i];
				value[i] = lower_bound_[i] + width * (bi + params[i]);
			}
			return true;
		}
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			cell_index = 0;
			for(size_t i = 0; i < DIM; ++i){
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				// Index cell_size_pref_sum = cell_sizes_.head(i).prod();
				Float bi = ( cell_index / cell_sizes_pref_sum_[i] ) % cell_sizes_[i];
				Float width = (upper_bound_[i]-lower_bound_[i])/(Float)cell_sizes_[i];
				params[i] = ( value[i] - lower_bound_[i] ) / width - bi;
				Float ci = (Index)params[i];
				cell_index += cell_sizes_pref_sum_[i] * ci;
				params[i] -= (Float)ci;
				// std::cout << i << " " << params[i] << " " << bi << " " << width << std::endl;
				assert( 0.0 < params[i] && params[i] < 1.0);
			}
			return true;
		}

		Float covering_radius(Index resl) const {
			Params width = (upper_bound_-lower_bound_) / cell_sizes_.template cast<Float>();
			return 0.5/(Float)(1<<resl) * sqrt((width*width).sum()); // norm2
		}
		Float neighbor_radius(Index resl) const {
			Params width = (upper_bound_-lower_bound_) / cell_sizes_.template cast<Float>();
			return 1.5/(Float)(1<<resl) * sqrt((width*width).sum()); // norm
		}

	 protected:
		~ScaleMap(){}
	};

	/// @brief Parameter Mapping Policy class which represents a discrete set of choices for 0 dimensional Nests
	/// @note NEST cell_size MUST agree with choices.size()
	template< int DIM, class Value, class Index, class Float >
	struct DiscreteChoiceMap {
		typedef Eigen::Array<Float,DIM,1> Params;		
		BOOST_STATIC_ASSERT_MSG(DIM==0,"DiscreteChoiceMap DIM must be == 0");
		std::vector<Value> choices;
		DiscreteChoiceMap(std::vector<Value> const & _choices) : choices(_choices){}
		///@brief sets value based only on cell_index
		///@note params has no meaning for zero-dimensional nests, only cell_index
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & /*params*/,
			Index cell_index,
			Value & value
		) const {
			if( cell_index >= choices.size() ) return false;
			value = choices[cell_index];
			return true;
		}
	 protected:
	 	~DiscreteChoiceMap(){}
	};


}
}

#endif
