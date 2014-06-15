#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include <Eigen/Dense>
#include <util/dilated_int.hh>
#include <iostream>
#include <vector>
#include <boost/static_assert.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

namespace scheme {
namespace nest {

	///@brief Base class for NEST
	///@tparam Index index type
	///@detail Base class for NEST, generic NEST interface
	template<class Index=size_t>
	struct NestBase {
		///@brief all nests must know their base_size
		NestBase(Index base_size) : base_size_(base_size) {}
		///@brief get the base_size of this NEST
		Index base_size() const { return base_size_; }
		///@brief generic virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@returns false if invalid index
		virtual	bool generic_set_state(
			std::vector<size_t> const & indices,
			Index base_index,
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
		size_t base_size_;
		void base_size(Index base_size) { base_size_ = base_size; }
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
	struct IdentityMap {
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Eigen::Matrix<Float,DIM,1> const & params,
			Value & value,
			Index base_index
		) const {
			for(size_t i = 0; i < DIM; ++i) value[i] = params[i];
				value[0] += base_index;
		}
	protected:
		~IdentityMap(){}
	};
	
	/// @brief Parameter Mapping Policy for cartesian grids
	/// @note NEST base_size MUST agree with base_sizes
	template< int DIM, class Value, class Index, class Float >
	struct ScaleMap {
		typedef Eigen::Matrix<Index,DIM,1> Indices;
		typedef Eigen::Matrix<Float,DIM,1> Params;		
		BOOST_STATIC_ASSERT_MSG(DIM>0,"ScaleMap DIM must be > 0");
		///@brief lower bound on value space		
		Params lower_bound;
		///@brief upper bound on value space base size 1
		Params upper_bound;
		///@brief distributes base_index accross dimensions
		Indices base_sizes;
		///@brief construct with default lb, ub, bs
		ScaleMap(){	base_sizes.fill(1);	lower_bound.fill(0); upper_bound.fill(1); }
		///@brief construct with default lb, ub
		ScaleMap(Indices const & bs) : base_sizes(bs) { lower_bound.fill(0); upper_bound.fill(1); }
		///@brief construct with default bs
		ScaleMap(Params const & lb, Params const & ub) : lower_bound(lb), upper_bound(ub) { base_sizes.fill(1); }
		///@brief construct with specified lb, ub and bs
		ScaleMap(Params const & lb, Params const & ub, Indices const & bs) : lower_bound(lb), upper_bound(ub), base_sizes(bs) {}
		///@brief sets value based on base_index and parameters using geometric bounds
		///@return false iff invalid parameters
		bool params_to_value(
			Eigen::Matrix<Float,DIM,1> const & params,
			Value & value,
			Index base_index
		) const {
			for(size_t i = 0; i < DIM; ++i){
				Float bi = ( base_index / base_sizes.head(i).prod() ) % base_sizes[i];
				value[i] = lower_bound[i] + (upper_bound[i]-lower_bound[i]) * (bi + params[i]);
			}
			return true;
		}
	protected:
		~ScaleMap(){}
	};

	/// @brief Parameter Mapping Policy class which represents a discrete set of choices for 0 dimensional Nests
	/// @note NEST base_size MUST agree with choices.size()
	template< int DIM, class Value, class Index, class Float >
	struct DiscreteChoiceMap {
		BOOST_STATIC_ASSERT_MSG(DIM==0,"DiscreteChoiceMap DIM must be == 0");
		std::vector<Value> choices;
		DiscreteChoiceMap(std::vector<Value> const & _choices) : choices(_choices){}
		///@brief sets value based only on base_index
		///@note params has no meaning for zero-dimensional nests, only base_index
		///@return false iff invalid parameters
		bool params_to_value(
			Eigen::Matrix<Float,DIM,1> const & /*params*/,
			Value & value,
			Index base_index
		) const {
			if( base_index >= choices.size() ) return false;
			value = choices[base_index];
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
		class Value = Eigen::Matrix<double,DIM,1>,
		template<int,class,class,class> class ParamMap = IdentityMap,
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
		NEST(Index base_size=1) : NestBase<Index>(base_size) {}
		///@brief for supporting ParamMaps, construct with bs
		NEST(Indices const & bs) : NestBase<Index>(bs.prod()), ParamMap<DIM,Value,Index,Float>(bs) {}
		///@brief for supporting ParamMaps, construct with default bs
		NEST(Params const & lb, Params const & ub) : ParamMap<DIM,Value,Index,Float>(lb,ub) {}
		///@brief for supporting ParamMaps, construct with specified lb, ub and bs
		NEST(Params const & lb, Params const & ub, Indices const & bs) : NestBase<Index>(bs.prod()), ParamMap<DIM,Value,Index,Float>(lb,ub,bs) {}

		///@brief get size of NEST at depth resl
		///@return size of NEST at resolution depth resl
		Index size(Index resl) const {
			return this->base_size() * ONE<<(DIM*resl);
		}

		///@brief set the state of this NEST to Value for index at depth resl
		///@return false iff invalid index
		bool set_state(Index index, Index resl){
			assert(index < size(resl));
			Index base_index = index >> (DIM*resl);
			Index hier_index = index & ((ONE<<(DIM*resl))-1);
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i){
				Index undilated = util::undilate<DIM>(hier_index>>i);
				params[i] = (static_cast<Float>(undilated) + 0.5 ) * scale;
			}
			return this->params_to_value( params, this->nonconst_value(), base_index );
		}

		///////////////////////////////////////
		//// generic interface functions
		///////////////////////////////////////

		///@brief generic virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@return false iff invalid index
		virtual bool generic_set_state(std::vector<size_t> const & hindices, Index base_index, size_t & iindex, Index resl) {
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i) params[i] = (static_cast<Float>(hindices[iindex]) + 0.5 ) * scale;
			iindex += DIM;
			return this->params_to_value( params, this->nonconst_value(), base_index );
		}
		///@brief virtual function returning size(resl)
		///@return size at depth resl
		virtual Index generic_size(Index resl) const { return size(resl); }
		///@brief virtual runction returning DIM
		///@return dimension of NEST
		virtual size_t generic_dim() const { return DIM; }
	};

	///@brief specialization of NEST for zero dimensional grids (discrete choices)
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
		NEST(Index base_size=1) : NestBase<Index>(base_size) {}
		NEST(std::vector<Value> const & _choices) : NestBase<Index>(_choices.size()), ParamMap<0,Value,Index,Float>(_choices) {}
		///@brief get num states at depth resl
		///@return number of status at depth resl
		Index size() const { return this->base_size(); }
		///@brief set state
		///@return false iff invalid index
		bool set_state(Index index){
			assert(index < this->base_size());
			Params params;
			return this->params_to_value( params, this->nonconst_value(), index );
		}

		///////////////////////////////////////
		//// generic interface functions
		///////////////////////////////////////

		///@brief generic virtual function to set the state of this nest
		///@detail will consume no indices from hindices vector, using only base_index
		///        for use in composite data structures containing NestBases
		virtual bool generic_set_state(std::vector<size_t> const & /*hindices*/, Index base_index, size_t & /*iindex*/, Index /*resl*/) {
			Params params;
			return this->params_to_value( params, this->nonconst_value(), base_index );
		}
		///@brief return size(resl)
		///@return base_size for these type
		virtual Index generic_size(Index /*resl*/) const { return this->base_size(); }
		///@brief get dimension of this nest
		///@return always 0 for these types
		virtual size_t generic_dim() const { return 0; }
	};

	void test(){ 
		using std::cout;
		using std::endl;

		Eigen::Vector2f lb(-1,-10),ub(1,2.4);
		Eigen::Matrix<size_t,2,1> bs(1,1);
		NEST<2,Eigen::RowVector2d,ScaleMap,StoreValue> nest(lb,ub,bs);
		size_t resl = 1;		
		cout << "test " << nest.size(resl) << endl;
		for(size_t i = 0; i < nest.size(resl); ++i){
			nest.set_state(i,resl);
			cout << "nest " << i << "i " << nest.value() << endl;
		}

		using namespace boost::assign;
		std::vector<double> choices; choices += 42.0,152.345,8049782.83402;
		NEST<0,double,DiscreteChoiceMap,StoreValue> nest0(choices);
		for(size_t i = 0; i < nest0.size(); ++i){
			nest0.set_state(i);
			std::cout << "nest0 " << i << "i " << nest0.value() << std::endl;
		}

	}


}
}

#endif
