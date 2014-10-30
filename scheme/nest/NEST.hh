#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/util/dilated_int.hh"
#include "scheme/util/StoragePolicy.hh"
#include "scheme/util/template_loop.hh"
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <boost/type_traits/make_signed.hpp>
#include <boost/any.hpp>

namespace scheme {
namespace nest {

	using util::StoreValue;
	using util::StorePointer;	

	///@brief Base class for NEST
	///@tparam Index index type
	///@detail Base class for NEST, virtual NEST interface
	template<class Index=size_t>
	struct NestBase {
		///@brief need virtual destructor
		virtual ~NestBase(){}
		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool 
		virtual_get_state(
			Index index, 
			Index resl, 
			boost::any & result
		) = 0;
		///@brief virtual virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@returns false if invalid index
		virtual	bool 
		virtual_get_state(
			std::vector<size_t> const & indices,
			Index cell_index,
			size_t & iindex,
			Index resl,
			boost::any & result
		) = 0;
		///@brief get the total size of this NEST at resl
		///@return number of possible states at depth resl
		virtual Index 
		virtual_size(Index resl) const = 0;
		///@brief get the dimension of this nest
		///@return dimension of Nest
		virtual size_t 
		virtual_dim() const = 0;
	};



	///////////////////////////////////////////////////////////////////////////////////////////////
	/// Main NEST template
	///////////////////////////////////////////////////////////////////////////////////////////////

	namespace maps { template< int DIM, class Value, class Index, class Float > struct UnitMap; }

	struct Empty {};

	///@brief templated nesting grids
	///@tparam DIM dimension of the grid
	///@tparam Value type of value the grid represents
	///@tparam ParamMap structure mapping from parameter space into Value space
	///@tparam StoragePolicy defines storage of Values, possibly allowing Nests to wrap pointers
	///@tparam Index index type, default size_t
	///@tparam Float floating point type for internal parameters, default double
	///@tparam bool is_virtual if you really want to optimize your code, you can set this to false
	///@note floats have plenty of precision for internal parameters
	template< int DIM,
		class Value = util::SimpleArray<DIM,double>,
		template<int,class,class,class> class ParamMap = maps::UnitMap,
		template<class> class StoragePolicy = StoreValue,
		class _Index = size_t,
		class Float = double,
		bool is_virtual = true
	>
	struct NEST : 
	    public boost::mpl::if_c<is_virtual,NestBase<_Index>,Empty>::type,
		public ParamMap<DIM,Value,_Index,Float>, 
	    public StoragePolicy<Value>
	{
		typedef _Index Index;
		typedef typename boost::make_signed<Index>::type SignedIndex;
		typedef NEST<DIM,Value,ParamMap,StoragePolicy,Index,Float,is_virtual> ThisType;
		typedef util::SimpleArray<DIM,Index> Indices;
		typedef util::SimpleArray<DIM,SignedIndex> SignedIndices;
		typedef util::SimpleArray<DIM,Float> Params;		
		typedef Value ValueType;
		typedef Index IndexType;
		typedef Float FloatType;
		typedef StoragePolicy<Value> StorageType;
		typedef ParamMap<DIM,Value,Index,Float> ParamMapType;

		static Index const ONE = 1;
		static int const DIMENSION = DIM;
		static Index const MAX_RESL_ONE_CELL = sizeof(Index)*8 / DIM;

		///@brief default ctor
		NEST() {}
		///@brief general constructor
		NEST(Index num_cells) : ParamMapType(num_cells) {}
		///@brief for supporting ParamMaps, construct with bs
		NEST(Indices const & bs) : ParamMapType(bs) {}
		///@brief for supporting ParamMaps, construct with default bs
		NEST(Params const & lb, Params const & ub) : ParamMapType(lb,ub) {}
		///@brief for supporting ParamMaps, construct with specified lb, ub and bs
		NEST(Params const & lb, Params const & ub, Indices const & bs) : ParamMapType(lb,ub,bs) {}

		///@brief get size of NEST at depth resl
		///@return size of NEST at resolution depth resl
		Index size(Index resl) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			return this->num_cells() * ONE<<(DIM*resl);
		}

		Index cell_index(Index index, Index resl) const {
			return index >> (DIM*resl);
		}

		///@brief set the input Value value for index at depth resl
		///@return false iff invalid index
		bool set_value(Index index, Index resl, Value & value) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			if(index >= size(resl)) return false;
			Index cell_index = index >> (DIM*resl);
			Index hier_index = index & ((ONE<<(DIM*resl))-1);
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i){
				Index undilated = util::undilate<DIM>(hier_index>>i);
				params[i] = (static_cast<Float>(undilated) + 0.5 ) * scale;
			}
			return this->params_to_value( params, cell_index, resl, value );
		}
		///@brief set the state of this NEST to Value for index at depth resl
		///@return false iff invalid index
		bool set_state(Index index, Index resl){
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			return set_value(index,resl,this->nonconst_value());
		}
		///@brief calls set_state(index,resl) then returns value()
		///@return value of set state
		Value const & set_and_get(Index index, Index resl){
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			bool set_state_returns_true = set_state(index,resl);
			assert( set_state_returns_true );
			return this->value();
		}
		///@brief return true iff index/resl is a valid bin
		bool check_state(Index index, Index resl) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			if(index >= size(resl)) return false;
			Value dummy;
			return set_value(index,resl,dummy);
		}
		///@brief get the index vector and cell index of the bin the Value is within
		///@returns true iff Value v is in a valid bin
		bool get_indicies(Value const & v, Index resl, Indices & indices_out, Index & cell_index_out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Params params;
			if( ! this->value_to_params( v, resl, params, cell_index_out ) ) return false;
			Float scale = Float(ONE<<resl);
			for(size_t i = 0; i < DIM; ++i){
				// this min for bounds check kinda sucks, but it's necessary if you want to allow points on the boundary
				indices_out[i] = std::min(((ONE<<resl)-1),static_cast<Index>(params[i]*scale));
			}
			return true;
		}
		///@brief get the index vector of a value WRT a particular cell, may be out of the cell bounds!
		///@detail this is used mainly for neighbor lookups -- some neighbors may be within the cell even if the value isn't
		///@returns nothing because the index vector isn't checked for validity
		void get_indicies_for_cell(Value const & v, Index resl, Index cell_index, Indices & indices_out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Params params;
			this->value_to_params_for_cell( v, resl, params, cell_index );
			Float scale = Float(ONE<<resl);
			for(size_t i = 0; i < DIM; ++i){
				// this crazy add/subtract avoids round towards 0 so params < 0 behave correctly
				Index const BIG = 12345678; //std::numeric_limits<Index>::max()/2;
				indices_out[i] = static_cast<Index>(params[i]*scale+BIG) - BIG;
			}
		}
		///@brief get the zorder index corresponding to and index vector and cell_index at resolution resl
		Index get_index(Indices const & indices, Index cell_index, Index resl) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i)	index |= util::dilate<DIM>(indices[i]) << i;
			index = index | (cell_index << (DIM*resl));
			return index;
		}
		///@brief get the zorder index for a Value v at resolution resl
		Index get_index(Value const & v, Index resl) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Index cell_index;
			Indices indices;
			if( !get_indicies(v,resl,indices,cell_index) ) return std::numeric_limits<Index>::max();
			return get_index(indices,cell_index,resl);
		}

		///@brief helper function for looping over neighbors and accumulating their indices
		template<class OutIter>
		void push_index(SignedIndices const & indices, Index cell_index, Index resl, OutIter out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i)	index |= util::dilate<DIM>(indices[i]) << i;
			index = index | (cell_index << (DIM*resl));
			*(out++) = index;
		}
		///@brief put the zorder indices of all neighbors of bin at (indices,cell_index) for resolution resl into OutIter out
		template<class OutIter>
		void get_neighbors(Indices const & indices, Index cell_index, Index resl, OutIter out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			// std::cout << indices.transpose() << std::endl;
			SignedIndices lb = ((indices.template cast<SignedIndex>()-(SignedIndex)1).max((SignedIndex)       0     ));
			SignedIndices ub = ((indices.template cast<SignedIndex>()+(SignedIndex)1).min((SignedIndex)((1<<resl)-1)));
			// std::cout << "IX " << indices.transpose() << " cell " << cell_index << std::endl;
			// std::cout << "LB " << lb.transpose() << std::endl;
			// std::cout << "UB " << ub.transpose() << std::endl;			
			boost::function<void(SignedIndices)> functor;
			functor = boost::bind( & ThisType::template push_index<OutIter>, this, _1, cell_index, resl, out );
			util::NESTED_FOR<DIM>(lb,ub,functor);
		}
		///@brief put the zorder indices of all neighbors of bin for Value v at resolution resl into OutIter out
		///@return false iff Value v itself dosen't have a valid index in this NEST
		template<class OutIter>
		bool get_neighbors(Value const & v, Index resl, OutIter out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			// get neighboring cells
			// get neighbors for each neighboring cell
			if( get_index(v,resl) == std::numeric_limits<Index>::max() ) return false;
			std::vector<Index> nbr_cells;
			std::back_insert_iterator<std::vector<Index> > inserter(nbr_cells);
			Float param_delta = 1.0 / (Float)(1<<resl);
			// std::cout << "DELTA " << param_delta << std::endl;
			this->get_neighboring_cells(v,resl,param_delta,inserter);
			BOOST_FOREACH( Index cell_index, nbr_cells ){
				// std::cout << "NBR_CELL " << cell_index << std::endl;
				Indices indices;
				this->get_indicies_for_cell(v,resl,cell_index,indices);
				this->get_neighbors(indices,cell_index,resl,out);
			}
			return true;
		}
		///@brief put the zorder indices of all neighbor bins within a particular cell for Value v at resolution
		///@brief resl into OutIter out
		template<class OutIter>
		void get_neighbors_for_cell(Value const & v, Index resl, Index cell_index, OutIter out) const {
			assert(resl<=MAX_RESL_ONE_CELL); // not rigerous check if Ncells > 1
			Indices indices;
			get_indicies_for_cell(v,resl,cell_index,indices);
			get_neighbors(indices,cell_index,resl,out);
		}

		///////////////////////////////////////
		//// virtual interface functions
		///////////////////////////////////////

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool virtual_get_state(Index index, Index resl, boost::any & result) {
			Value & v = *boost::any_cast<Value*>(result);
			bool status = set_value( index, resl, v );
			if(status) this->set_stored_value(v);
			return status;
		}
		///@brief virtual virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@return false iff invalid index
		virtual bool virtual_get_state(
			std::vector<size_t> const & hindices,
			Index cell_index,
			size_t & iindex,
			Index resl,
			boost::any & result
		) {
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i) params[i] = (static_cast<Float>(hindices[iindex+i]) + 0.5 ) * scale;
			iindex += DIM;
			Value & v( *boost::any_cast<Value*>(result) );
			bool status = this->params_to_value( params, cell_index, resl, v );
			if(status) this->set_stored_value(v);
			return status;
		}
		///@brief virtual function returning size(resl)
		///@return size at depth resl
		virtual Index virtual_size(Index resl) const { return size(resl); }
		///@brief virtual runction returning DIM
		///@return dimension of NEST
		virtual size_t virtual_dim() const { return DIM; }
	};

	////////////////////////////////////////////////////////////////////////////////
	///@brief specialization of NEST for zero dimensional grids (discrete choices)
	///@detail some logic has to be shortcurcited to work with 0-dimensional grids
	////////////////////////////////////////////////////////////////////////////////
	template<
		class Value,
		template<int,class,class,class> class ParamMap,
		template<class> class StoragePolicy,
		class _Index,
		class Float,
		bool is_virtual
	>
	struct NEST<0,Value,ParamMap,StoragePolicy,_Index,Float,is_virtual> : 
	    		  public boost::mpl::if_c<is_virtual,NestBase<_Index>,Empty>::type,
	              public ParamMap<0,Value,_Index,Float>, 
	              public StoragePolicy<Value>
	{
		typedef _Index Index;
		typedef char Params; // no params
		///@brief choices vector constructor
		NEST(std::vector<Value> const & _choices) : ParamMap<0,Value,Index,Float>(_choices){}
		///@brief get num states at depth resl
		///@return number of status at depth resl
		Index size(Index /*resl*/=0) const { return this->num_cells(); }
		///@brief set state
		///@return false iff invalid index
		bool set_state(Index index, Index resl, Value & val){
			// assert(index < this->num_cells());
			Params params;
			return this->params_to_value( params, index, resl, val );
		}
		///@brief set state
		///@return false iff invalid index
		bool set_state(Index index, Index resl=0){
			// assert(index < this->num_cells());
			Params params;
			return this->params_to_value( params, index, resl, this->nonconst_value() );
		}
		///@brief calls set_state(index,resl) then returns value()
		///@return value of set state
		Value const & set_and_get(Index index, Index resl=0){
			assert( set_state(index,resl) );
			return this->value();
		}

		///////////////////////////////////////
		//// virtual interface functions
		///////////////////////////////////////

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool virtual_get_state(Index index, Index resl, boost::any & result) {
			Value & v( *boost::any_cast<Value*>(result) );
			bool status = set_state( index, resl, v );
			if(status) this->set_stored_value(v);
			return status;
		}

		///@brief virtual virtual function to set the state of this nest
		///@detail will consume no indices from hindices vector, using only cell_index
		///        for use in composite data structures containing NestBases
		virtual
		bool
		virtual_get_state(
			std::vector<size_t> const & /*hindices*/,
			Index cell_index,
			size_t & /*iindex*/,
			Index resl,
			boost::any & result
		) {
			Params params;
			Value & v( *boost::any_cast<Value*>(result) );
			bool status = this->params_to_value( params, cell_index, resl, v );
			if(status) this->set_stored_value(v);
			return status;
		}
		///@brief return size(resl)
		///@return num_cells for these type
		virtual Index virtual_size(Index /*resl*/) const { return this->num_cells(); }
		///@brief get dimension of this nest
		///@return always 0 for these types
		virtual size_t virtual_dim() const { return 0; }
	};


}
}

#endif
