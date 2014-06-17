#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include <Eigen/Dense>
#include <util/dilated_int.hh>
#include <iostream>
#include <vector>
#include <boost/static_assert.hpp>

namespace scheme {
namespace nest {

	using namespace Eigen;
	typedef Eigen::Matrix<size_t,1,1> Vector1s;
	typedef Eigen::Matrix<size_t,2,1> Vector2s;
	typedef Eigen::Matrix<size_t,3,1> Vector3s;
	typedef Eigen::Matrix<size_t,4,1> Vector4s;			
	typedef Eigen::Matrix<size_t,5,1> Vector5s;			
	typedef Eigen::Matrix<size_t,6,1> Vector6s;			

	///@brief Base class for NEST
	///@tparam Index index type
	///@detail Base class for NEST, generic NEST interface
	template<class Index=size_t>
	struct NestBase {
		///@brief all nests must know their cell_size
		NestBase(Index cell_size) : cell_size_(cell_size) {}
		///@brief need virtual destructor
		virtual ~NestBase(){}
		///@brief get the cell_size of this NEST
		Index cell_size() const { return cell_size_; }
		///@brief generic virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@returns false if invalid index
		virtual	bool generic_set_state(
			std::vector<size_t> const & indices,
			Index cell_index,
			size_t & iindex,
			Index resl
		) = 0;
		///@brief get the total size of this NEST at resl
		///@return number of possible states at depth resl
		virtual Index generic_size(Index resl) const = 0;
		///@brief get the dimension of this nest
		///@return dimension of Nest
		virtual size_t generic_dim() const = 0;
	protected:
		size_t cell_size_;
		void cell_size(Index cell_size) { cell_size_ = cell_size; }
	};



	/// @brief Storage Policy Class, store by value
	/// @tparam Value ValueType stored
	template< class Value >
	struct StoreValue {
		///@return const reference to stored Value
		Value const & value() const { return value_; }
	protected:
		///@return nonconst reference to stored Value
		Value & nonconst_value() { return value_; }
		Value value_;
		~StoreValue(){}
	};

	/// @brief Store-by-pointer policy
	/// @tparam Value ValueType stored
	/// @note addes the ability to set the pointer
	template< class Value >
	struct StorePointer {
		///@return const reference to stored Value
		Value const & value() const { return *value_; }
		/// @brief switch the pointer the this policy manages
		/// @param new_pointer 
		void set_pointer(Value * new_pointer) { value_ = new_pointer; }
	protected:
		///@return nonconst reference to stored Value
		Value & nonconst_value() { return *value_; }
		Value * value_;
		~StorePointer(){}
	};

	/// @brief Parameter to Value Map Policy Class
	template< int DIM, class Value, class Index, class Float >
	struct UnitMap {
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Eigen::Matrix<Float,DIM,1> const & params,
			Index cell_index,
			Value & value
		) const {
			for(size_t i = 0; i < DIM; ++i) value[i] = params[i];
			value[0] += (Float)cell_index;
			return true;
		}
		bool value_to_params(
			Value const & value,
			Eigen::Matrix<Float,DIM,1> & params,
			Index & cell_index
		) const {
			for(size_t i = 0; i < DIM; ++i) params[i] = value[i];
			cell_index = (Index)value[0];
			params[0] -= (Float)cell_index;
			return true;
		}
	protected:
		~UnitMap(){}
	};
	
	/// @brief Parameter Mapping Policy for cartesian grids
	/// @note NEST cell_size MUST agree with cell_sizes_
	template< int DIM, class Value, class Index, class Float >
	struct ScaleMap {
		typedef Eigen::Matrix<Index,DIM,1> Indices;
		typedef Eigen::Matrix<Float,DIM,1> Params;		
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
			Eigen::Matrix<Float,DIM,1> const & params,
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
			Eigen::Matrix<Float,DIM,1> & params,
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

	protected:
		~ScaleMap(){}
	};

	/// @brief Parameter Mapping Policy class which represents a discrete set of choices for 0 dimensional Nests
	/// @note NEST cell_size MUST agree with choices.size()
	template< int DIM, class Value, class Index, class Float >
	struct DiscreteChoiceMap {
		BOOST_STATIC_ASSERT_MSG(DIM==0,"DiscreteChoiceMap DIM must be == 0");
		std::vector<Value> choices;
		DiscreteChoiceMap(std::vector<Value> const & _choices) : choices(_choices){}
		///@brief sets value based only on cell_index
		///@note params has no meaning for zero-dimensional nests, only cell_index
		///@return false iff invalid parameters
		bool params_to_value(
			Eigen::Matrix<Float,DIM,1> const & /*params*/,
			Index cell_index,
			Value & value
		) const {
			if( cell_index >= choices.size() ) return false;
			value = choices[cell_index];
			return true;
		}
	};

	///@brief templated nesting grids
	///@tparam DIM dimension of the grid
	///@tparam Value type of value the grid represents
	///@tparam ParamMap structure mapping from parameter space into Value space
	///@tparam StoragePolicy defines storage of Values, possibly allowing Nests to wrap pointers
	///@tparam Index index type, default size_t
	///@tparam Float floating point type for internal parameters, default float
	///@note floats have plenty of precision for internal parameters
	template< int DIM,
		class Value = Eigen::Matrix<float,DIM,1>,
		template<int,class,class,class> class ParamMap = UnitMap,
		template<class> class StoragePolicy = StoreValue,
		class Index = size_t,
		class Float = float
	>
	struct NEST : 
	    public NestBase<Index>,
		public ParamMap<DIM,Value,Index,Float>, 
	    public StoragePolicy<Value>
	{
		typedef Eigen::Matrix<Index,DIM,1> Indices;
		typedef Eigen::Matrix<Float,DIM,1> Params;		
		static Index const ONE = 1;

		///@brief default ctor
		NEST() : NestBase<Index>(1) {}
		///@brief general constructor
		NEST(Index cell_size) : NestBase<Index>(cell_size) {}
		///@brief for supporting ParamMaps, construct with bs
		NEST(Indices const & bs) : NestBase<Index>(bs.prod()), ParamMap<DIM,Value,Index,Float>(bs) {}
		///@brief for supporting ParamMaps, construct with default bs
		NEST(Params const & lb, Params const & ub) : ParamMap<DIM,Value,Index,Float>(lb,ub) {}
		///@brief for supporting ParamMaps, construct with specified lb, ub and bs
		NEST(Params const & lb, Params const & ub, Indices const & bs) : 
			NestBase<Index>(bs.prod()), ParamMap<DIM,Value,Index,Float>(lb,ub,bs) {}

		///@brief get size of NEST at depth resl
		///@return size of NEST at resolution depth resl
		Index size(Index resl) const {
			return this->cell_size() * ONE<<(DIM*resl);
		}

		///@brief set the state of this NEST to Value for index at depth resl
		///@return false iff invalid index
		bool set_state(Index index, Index resl){
			if(index >= size(resl)) return false;
			Index cell_index = index >> (DIM*resl);
			Index hier_index = index & ((ONE<<(DIM*resl))-1);
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i){
				Index undilated = util::undilate<DIM>(hier_index>>i);
				params[i] = (static_cast<Float>(undilated) + 0.5 ) * scale;
			}
			return this->params_to_value( params, cell_index, this->nonconst_value() );
		}

		///@brief calls set_state(index,resl) then returns value()
		///@return value of set state
		Value const & set_and_get(Index index, Index resl){
			assert( set_state(index,resl) );
			return this->value();
		}

		Index get_index(Value const & v, Index resl) const{
			Index cell_index;
			Params params;
			if( ! this->value_to_params( v, params, cell_index ) ) return std::numeric_limits<Index>::max();
			Indices indices;
			Float scale = Float(ONE<<resl);
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i){
				Index undilated = static_cast<Index>(params[i]*scale);
				index |= util::dilate<DIM>(undilated) << i;
			}
			// std::cout << v[0] <<"v "<< resl <<"r " << cell_index << "bi " << index << "i " << (cell_index << (DIM*resl)) <<  std::endl;
			index = index | (cell_index << (DIM*resl));
			return index;

		}

		///////////////////////////////////////
		//// generic interface functions
		///////////////////////////////////////

		///@brief generic virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@return false iff invalid index
		virtual bool generic_set_state(std::vector<size_t> const & hindices, Index cell_index, size_t & iindex, Index resl) {
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i) params[i] = (static_cast<Float>(hindices[iindex]) + 0.5 ) * scale;
			iindex += DIM;
			return this->params_to_value( params, cell_index, this->nonconst_value() );
		}
		///@brief virtual function returning size(resl)
		///@return size at depth resl
		virtual Index generic_size(Index resl) const { return size(resl); }
		///@brief virtual runction returning DIM
		///@return dimension of NEST
		virtual size_t generic_dim() const { return DIM; }
	};

	////////////////////////////////////////////////////////////////////////////////
	///@brief specialization of NEST for zero dimensional grids (discrete choices)
	///@detail some logic has to be shortcurcited to work with 0-dimensional grids
	////////////////////////////////////////////////////////////////////////////////
	template<
		class Value,
		template<int,class,class,class> class ParamMap,
		template<class> class StoragePolicy,
		class Index,
		class Float
	>
	struct NEST<0,Value,ParamMap,StoragePolicy,Index,Float> : 
	              public NestBase<Index>,
	              public ParamMap<0,Value,Index,Float>, 
	              public StoragePolicy<Value>
	{
		typedef Eigen::Matrix<Float,0,1> Params;		
		NEST(Index cell_size=1) : NestBase<Index>(cell_size) {}
		NEST(std::vector<Value> const & _choices) : NestBase<Index>(_choices.size()), ParamMap<0,Value,Index,Float>(_choices) {}
		///@brief get num states at depth resl
		///@return number of status at depth resl
		Index size(Index /*resl*/=0) const { return this->cell_size(); }
		///@brief set state
		///@return false iff invalid index
		bool set_state(Index index, Index /*resl*/=0){
			// assert(index < this->cell_size());
			Params params;
			return this->params_to_value( params, index, this->nonconst_value() );
		}

		///@brief calls set_state(index,resl) then returns value()
		///@return value of set state
		Value const & set_and_get(Index index, Index resl=0){
			assert( set_state(index,resl) );
			return this->value();
		}

		///////////////////////////////////////
		//// generic interface functions
		///////////////////////////////////////

		///@brief generic virtual function to set the state of this nest
		///@detail will consume no indices from hindices vector, using only cell_index
		///        for use in composite data structures containing NestBases
		virtual bool generic_set_state(std::vector<size_t> const & /*hindices*/, Index cell_index, size_t & /*iindex*/, Index /*resl*/) {
			Params params;
			return this->params_to_value( params, cell_index, this->nonconst_value() );
		}
		///@brief return size(resl)
		///@return cell_size for these type
		virtual Index generic_size(Index /*resl*/) const { return this->cell_size(); }
		///@brief get dimension of this nest
		///@return always 0 for these types
		virtual size_t generic_dim() const { return 0; }
	};


}
}

#endif
