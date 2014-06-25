#ifndef INCLUDED_scheme_nest_maps_parameter_maps_HH
#define INCLUDED_scheme_nest_maps_parameter_maps_HH

#include <util/template_loop.hh>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace maps {

	using namespace Eigen;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///	Parameter Mappings
	///////////////////////////////////////////////////////////////////////////////////////////////

	///@brief Parameter to Value Map Policy Class Concept
	///@detail JUST A DEFINITION OF THE CONCEPT AND PLACEHOLDER!
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
	struct ParamMap {
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
		) const ;

		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const ;

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
	struct UnitMap : public ParamMap<DIM,Value,Index,Float> {
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
			value_to_params_for_cell(value,params,0);
			cell_index = (Index)value[0];
			params[0] -= (Float)cell_index;
			return true;
		}
		///@brief get params repr of Value wrt cell cell_index
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Params & params,
			Index cell_index
		) const {
			for(size_t i = 0; i < DIM; ++i) params[i] = value[i];
			params[0] -= (Float)cell_index;
		}
		///@brief return the cell_index of neighboring cells within delta of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float param_delta, OutIter out) const {
			// this BIG thing is to ensure rounding goes down
			int const BIG = 12345678;
			int lb = std::max(                0, static_cast<int>( value[0]-param_delta + BIG ) - BIG );
			int ub = std::min((int)cell_size_-1, static_cast<int>( value[0]+param_delta + BIG ) - BIG );
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
	///@note bounds and cell indices are represented as Eigen Arrays (like params) NOT Value Types
	template<
		int DIM,
		class Value=Eigen::Matrix<double,DIM,1>,
		class Index=size_t,
		class Float=double
	>
	struct ScaleMap : ParamMap<DIM,Value,Index,Float> {
		static int const DIMENSION = DIM;
		typedef ScaleMap<DIM,Value,Index,Float> ThisType;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<int,DIM,1> SignedIndices;
		typedef Eigen::Array<Float,DIM,1> Params;		
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
	 private:
		///@brief lower bound on value space		
		Params lower_bound_;
		///@brief upper bound on value space base size 1
		Params upper_bound_;
		///@brief distributes cell_index accross dimensions
		Params cell_width_;
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
			for(size_t i = 0; i < DIM; ++i){
				cell_sizes_pref_sum_[i] = cell_sizes_.head(i).prod();
				cell_width_[i] = (upper_bound_[i]-lower_bound_[i])/(Float)cell_sizes_[i];
				assert(upper_bound_[i] > lower_bound_[i]);
			}
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
				value[i] = lower_bound_[i] + cell_width_[i] * (bi + params[i]);
			}
			return true;
		}

		///@brief sets params/cell_index from value
		bool value_to_params(
			Value const & value,
			Params & params,
			Index & cell_index
		) const {
			value_to_params_for_cell(value,params,0);
			cell_index = 0;
			for(size_t i = 0; i < DIM; ++i){
				assert( params[i] >= 0.0 );
				assert( params[i] < cell_sizes_[i]+1);
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				// Index cell_size_pref_sum = cell_sizes_.head(i).prod();
				Float ci = (Index)params[i];
				cell_index += cell_sizes_pref_sum_[i] * ci;
				params[i] -= (Float)ci;
				assert( 0.0 <= params[i] && params[i] < 1.0);
			}
			return true;
		}

		///@brief sets params/cell_index from value
		void value_to_params_for_cell(
			Value const & value,
			Params & params,
			Index cell_index
		) const {
			Indices cell_indices = cellindex_to_indices(cell_index);
			// std::cout << "CIDX " << cell_indices.transpose() << std::endl;
			for(size_t i = 0; i < DIM; ++i){
				assert(cell_sizes_[i] > 0);
				assert(cell_sizes_[i] < 100000);
				assert(lower_bound_[i] < upper_bound_[i]);
				// Index cell_size_pref_sum = cell_sizes_.head(i).prod();
				params[i] = ( value[i] - lower_bound_[i] ) / cell_width_[i];
				params[i] -= (Float)cell_indices[i];
			}
		}

		Index indices_to_cellindex(Indices const & indices) const {
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i){
				assert(indices[i] < cell_sizes_[i]);
				index += indices[i]*cell_sizes_pref_sum_[i];
			}
			return index;
		}

		Indices cellindex_to_indices(Index index) const {
			Indices indices;						
			for(size_t i = 0; i < DIM; ++i){
				indices[i] = ( index / cell_sizes_pref_sum_[i] ) % cell_sizes_[i];
				assert(indices[i] < cell_sizes_[i]);
			}
			return indices;
		}

		template<class OutIter>
		void push_cell_index(SignedIndices const & indices, OutIter out) const {
			*(out++) = indices_to_cellindex(indices.template cast<size_t>());
		}

		// template<class OutIter>
		// void get_neighbors(Indices const & indices, Index cell_index, Index resl, OutIter out)  {
		// 	// std::cout << indices.transpose() << std::endl;
		// 	SignedIndices lb = ((indices.template cast<int>()-1).max(     0     ));
		// 	SignedIndices ub = ((indices.template cast<int>()+1).min((1<<resl)-1));
		// 	// std::cout << "IX " << indices.transpose() << " cell " << cell_index << std::endl;
		// 	// std::cout << "LB " << lb.transpose() << std::endl;
		// 	// std::cout << "UB " << ub.transpose() << std::endl;			
		// 	boost::function<void(SignedIndices)> functor;
		// 	functor = boost::bind( & ThisType::template push_index<OutIter>, this, _1, cell_index, resl, out );
		// 	util::NESTED_FOR<DIM>(lb,ub,functor);
		// }

		///@brief return the cell_index of neighboring cells within delta of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float param_delta, OutIter out) const {
			assert( param_delta > 0);
			// convert to value space, decided against this
			// Params delta_param;
			// for(size_t i = 0; i < DIM; ++i) delta_param[i] = delta / cell_width_[i];
			Params params;
			value_to_params_for_cell(value,params,0);
			int const BIG = 12345678;
			SignedIndices lb = (params-param_delta+(Float)BIG).max((Float)BIG).template cast<int>() - BIG ;
			SignedIndices ub = (params+param_delta).template cast<Index>().min(cell_sizes_-1).template cast<int>();
			// std::cout << "PM " << params.transpose() << std::endl;
			// std::cout << "DL " << delta_param.transpose() << std::endl;
			// std::cout << "LB " << lb.transpose() << std::endl;
			// std::cout << "UB " << ub.transpose() << std::endl;
			boost::function<void(SignedIndices)> functor;
			functor = boost::bind( & ThisType::template push_cell_index<OutIter>, this, _1, out );
			util::NESTED_FOR<DIM>(lb,ub,functor);
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
}

#endif
