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
	///@tparam DIM the dimension number of the input parameter space
	///@tparam Value the output value type, default Eigen Matrix
	///@tparam Index index type, default size_t
	///@tparam Float float type, default double
	template<
		int DIM,
		class Value=Eigen::Matrix<double,DIM,1>,
		class Index=size_t,
		class Float=double
	>
	struct UnitMap {
		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<Float,DIM,1> Params;

		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
		///@brief constructor
		UnitMap(Index cell_size=1) : cell_size_(cell_size) {}
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
		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			value_to_params_unitcell(value,params);
			cell_index = (Index)value[0];
			params[0] -= (Float)cell_index;
			return true;
		}
		///@brief for unit cell
		///@note necessary only for neighbor lookup		
		void value_to_params_unitcell(
			Value const & value,
			Params & params
		) const {
			for(size_t i = 0; i < DIM; ++i) params[i] = value[i];
		}
		///@brief return the cell_index of neighboring cells within radius of value
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float radius, OutIter out) {
			// this BIG thing is to ensure rounding goes down
			Index const BIG = 12345678;
			int lb = static_cast<int>( fmax(     0.0      ,value[0]-radius) + BIG ) - BIG;
			int ub = static_cast<int>( fmin(cell_size_+0.5,value[0]+radius) + BIG ) - BIG;
			// std::cout << "lb " << lb << " ub "  << ub << std::endl;
			// assert(lb<=ub);
			for(int i = lb; i <= ub; ++i) *(out++) = i;
		}
		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const { return 0.5/(Float)(1<<resl) * sqrt(DIM); }
		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const { return 1.5/(Float)(1<<resl); }
		///@brief cell size
		Index cell_size() const { return cell_size_; }
		virtual ~UnitMap(){}
	 private:
	 	///@brief number of cells
		Index cell_size_;
	};
	


	///@brief Parameter Mapping Policy for cartesian grids
	///@tparam DIM the dimension number of the input parameter space
	///@tparam Value the output value type, default Eigen Matrix
	///@tparam Index index type, default size_t
	///@tparam Float float type, default double
	///@note NEST cell_size MUST agree with cell_sizes_
	template<
		int DIM,
		class Value=Eigen::Matrix<double,DIM,1>,
		class Index=size_t,
		class Float=double
	>
	struct ScaleMap {
		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
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
		Index cell_size_;
	 public:
		///@brief construct with default lb, ub, bs
		ScaleMap(){	cell_sizes_.fill(1); lower_bound_.fill(0); upper_bound_.fill(1); init(); }
		///@brief construct with default lb, ub
		ScaleMap(Indices const & bs) : 
			cell_sizes_(bs) { lower_bound_.fill(0); upper_bound_.fill(1); init(); }
		///@brief construct with default bs
		ScaleMap(Params const & lb, Params const & ub) : 
			lower_bound_(lb), upper_bound_(ub) { cell_sizes_.fill(1); init(); }
		///@brief construct with specified lb, ub and bs
		ScaleMap(Params const & lb, Params const & ub, Indices const & bs) : 
			lower_bound_(lb), upper_bound_(ub), cell_sizes_(bs) { init(); }

		///@brief sets up cell_size_pref_sum
		void init(){
			cell_size_ = cell_sizes_.prod();
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

		///@brief sets params/cell_index from value
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			value_to_params_unitcell(value,params);
			cell_index = 0;
			for(size_t i = 0; i < DIM; ++i){
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				// Index cell_size_pref_sum = cell_sizes_.head(i).prod();
				Float ci = (Index)params[i];
				cell_index += cell_sizes_pref_sum_[i] * ci;
				params[i] -= (Float)ci;
				// std::cout << i << " " << params[i] << " " << bi << " " << width << std::endl;
				assert( 0.0 < params[i] && params[i] < 1.0);
			}
			return true;
		}

		///@brief sets params/cell_index from value
		void value_to_params_unitcell(
			Value const & value,
			Params & params
		) const {
			for(size_t i = 0; i < DIM; ++i){
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				// Index cell_size_pref_sum = cell_sizes_.head(i).prod();
				Float width = (upper_bound_[i]-lower_bound_[i])/(Float)cell_sizes_[i];
				params[i] = ( value[i] - lower_bound_[i] ) / width;
			}
		}

		///@brief return the cell_index of neighboring cells within radius of value
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float radius, OutIter out) {
			// this BIG thing is to ensure rounding goes down
			Index const BIG = 12345678;
		}

		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const {
			Params width = (upper_bound_-lower_bound_) / cell_sizes_.template cast<Float>();
			return 0.5/(Float)(1<<resl) * sqrt((width*width).sum()); // norm2
		}

		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const {
			Params width = (upper_bound_-lower_bound_) / cell_sizes_.template cast<Float>();
			return 1.5/(Float)(1<<resl) * width.minCoeff(); // norm
		}

		///@brief cell size
		Index cell_size() const { return cell_size_; }
	 	virtual ~ScaleMap(){}
	};

	/// @brief Parameter Mapping Policy class which represents a discrete set of choices for 0 dimensional Nests
	/// @note NEST cell_size MUST agree with choices.size()
	template< int DIM, class Value, class Index, class Float >
	struct DiscreteChoiceMap {
		BOOST_STATIC_ASSERT_MSG(DIM==0,"DiscreteChoiceMap DIM must be == 0");
		static int const DIMENSION = DIM;
		typedef Eigen::Array<Float,DIM,1> Params;		
		std::vector<Value> choices;
		DiscreteChoiceMap(std::vector<Value> const & _choices) : choices(_choices), cell_size_(choices.size()) {}
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
		Index cell_size() const { return cell_size_; }
		virtual ~DiscreteChoiceMap(){}
	 private:
	 	Index cell_size_;
	};


}
}

#endif
