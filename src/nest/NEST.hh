#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include <Eigen/Dense>
#include <util/dilated_int.hh>
#include <iostream>

namespace scheme {
namespace nest {

	template< class Value >
	struct StoreValue {
		Value const & value() const { return value_; }
	protected:
		Value & value() { return value_; }
		Value value_;
		~StoreValue(){}
	};
	template< class Value >

	struct StorePointer {
		Value const & value() const { return *value_; }
		void set_pointer(Value * new_pointer) { value_ = new_pointer; }
	protected:
		Value & value() { return *value_; }
		Value * value_;
		~StorePointer(){}
	};

	template< int DIM, class Value, class Float >
	struct IdentityMap {
		void params_to_value(Eigen::Array<Float,DIM,1> const & params, Value & value) const {
			for(size_t i = 0; i < DIM; ++i) value[i] = params[i];
		}
	protected:
		~IdentityMap(){}
	};
	
	template< int DIM, class Value, class Float >
	struct ScaleMap {
		Eigen::Array<Float,DIM,1> lower_bound;
		Eigen::Array<Float,DIM,1> upper_bound;
		void params_to_value(Eigen::Array<Float,DIM,1> const & params, Value & value) const {
			for(size_t i = 0; i < DIM; ++i){
				value[i] = lower_bound[i]+ (upper_bound[i]-lower_bound[i])*params[i];
			}
		}
	protected:
		~ScaleMap(){}
	}; 

	template< int DIM,
		class Float   = double, 
		class Value   = Eigen::Array<Float,DIM,1>,
		template<int,class,class> class ParamMap = IdentityMap,
		template<class> class StoragePolicy = StoreValue,
		class Index   = uint64_t
	>
	struct NEST : public ParamMap<DIM,Value,Float>, 
	              public StoragePolicy<Value>
	{
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<Float,DIM,1> Params;		
		static Index const ONE = 1;
		Eigen::Array<Index,DIM,1> base_size;
		NEST(){
			base_size.fill(1);
		}
		Index size(Index resl) const {
			return base_size.prod() * ONE<<(DIM*resl);
		}
		Indices expand_index(Index index, Index resl){
			assert(index < size(resl));
			Index base_index = index >> (Index)DIM*resl;
			Index hier_index = index & ((ONE<<((Index)3*resl))-1);
			Indices indices;
			for(size_t i = 0; i < DIM; ++i){
				Index bi = ( base_index / base_size.head(i).prod() ) % base_size[i];
				Index hi = util::undilate<DIM>(hier_index>>i);
				indices[i] = (bi<<resl) + hi;
			}
			return indices;
		}
		Value const & set_state(Index index, Index resl){
			Params params = expand_index(index,resl).template cast<Float>();
			params = (params + 0.5) / Float(ONE<<resl);
			this->params_to_value( params, this->value() );
			return this->value();
		}
	};

	void test(){ 
		using std::cout;
		using std::endl;
		NEST<2,float,Eigen::Vector2f,ScaleMap,StoreValue> nest;
		nest.base_size << 1,1;
		nest.lower_bound << -1,-1;
		nest.upper_bound <<  3, 2;		
		cout << "test " << nest.size(4) << " / " << nest.base_size.transpose() << endl;
		for(size_t i = 0; i < nest.size(4); ++i){
			cout << "INDEX " << nest.set_state(i,4).transpose() << endl;
		}
	}


}
}

#endif
