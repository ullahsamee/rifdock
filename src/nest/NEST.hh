#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include <Eigen/Dense>
#include <util/dilated_int.hh>
#include <vector>
#include <util/storage_policy.hh>
#include <util/template_loop.hh>
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace scheme {
namespace nest {

	
	using util::StoreValue;
	using util::StorePointer;	

	using namespace Eigen;
	typedef Eigen::Matrix<size_t,1,1> Vector1s;
	typedef Eigen::Matrix<size_t,2,1> Vector2s;
	typedef Eigen::Matrix<size_t,3,1> Vector3s;
	typedef Eigen::Matrix<size_t,4,1> Vector4s;			
	typedef Eigen::Matrix<size_t,5,1> Vector5s;			
	typedef Eigen::Matrix<size_t,6,1> Vector6s;			

	///@brief Base class for NEST
	///@tparam Index index type
	///@detail Base class for NEST, virtual NEST interface
	template<class Index=size_t>
	struct NestBase {
		///@brief all nests must know their cell_size
		NestBase(Index cell_size=1) : cell_size_(cell_size) {}
		///@brief need virtual destructor
		virtual ~NestBase(){}
		///@brief get the cell_size of this NEST
		Index cell_size() const { return cell_size_; }

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool virtual_set_state(Index index, Index resl) = 0;

		///@brief virtual virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@returns false if invalid index
		virtual	bool virtual_set_state(
			std::vector<size_t> const & indices,
			Index cell_index,
			size_t & iindex,
			Index resl
		) = 0;
		///@brief get the total size of this NEST at resl
		///@return number of possible states at depth resl
		virtual Index virtual_size(Index resl) const = 0;
		///@brief get the dimension of this nest
		///@return dimension of Nest
		virtual size_t virtual_dim() const = 0;
	protected:
		size_t cell_size_;
		void cell_size(Index cell_size) { cell_size_ = cell_size; }
	};



	///////////////////////////////////////////////////////////////////////////////////////////////
	/// Main NEST template
	///////////////////////////////////////////////////////////////////////////////////////////////

	template< int DIM, class Value, class Index, class Float > struct UnitMap;

	///@brief templated nesting grids
	///@tparam DIM dimension of the grid
	///@tparam Value type of value the grid represents
	///@tparam ParamMap structure mapping from parameter space into Value space
	///@tparam StoragePolicy defines storage of Values, possibly allowing Nests to wrap pointers
	///@tparam Index index type, default size_t
	///@tparam Float floating point type for internal parameters, default double
	///@note floats have plenty of precision for internal parameters
	template< int DIM,
		class Value = Eigen::Matrix<double,DIM,1>,
		template<int,class,class,class> class ParamMap = UnitMap,
		template<class> class StoragePolicy = StoreValue,
		class Index = size_t,
		class Float = double
	>
	struct NEST : 
	    public NestBase<Index>,
		public ParamMap<DIM,Value,Index,Float>, 
	    public StoragePolicy<Value>
	{
		typedef NEST<DIM,Value,ParamMap,StoragePolicy,Index,Float> ThisType;
		typedef Eigen::Array<Index,DIM,1> Indices;
		typedef Eigen::Array<int,DIM,1> SignedIndices;
		typedef Eigen::Array<Float,DIM,1> Params;		
		typedef Value ValueType;
		typedef Index IndexType;
		typedef Float FloatType;
		typedef StoragePolicy<Value> StorageType;
		typedef ParamMap<DIM,Value,Index,Float> ParamMapType;

		static Index const ONE = 1;
		static int const DIMENSION = DIM;

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

		///@brief set the input Value value for index at depth resl
		///@return false iff invalid index
		bool set_value(Index index, Index resl, Value & value) const {
			if(index >= size(resl)) return false;
			Index cell_index = index >> (DIM*resl);
			Index hier_index = index & ((ONE<<(DIM*resl))-1);
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i){
				Index undilated = util::undilate<DIM>(hier_index>>i);
				params[i] = (static_cast<Float>(undilated) + 0.5 ) * scale;
			}
			return this->params_to_value( params, cell_index, value );
		}
		///@brief set the state of this NEST to Value for index at depth resl
		///@return false iff invalid index
		bool set_state(Index index, Index resl){
			return set_value(index,resl,this->nonconst_value());
		}
		///@brief calls set_state(index,resl) then returns value()
		///@return value of set state
		Value const & set_and_get(Index index, Index resl){
			bool tf = set_state(index,resl);
			assert( tf );
			return this->value();
		}

		bool check_state(Index index, Index resl) const {
			if(index >= size(resl)) return false;
			Value dummy;
			return set_value(index,resl,dummy);
		}

		bool get_indicies(Value const & v, Index resl, Indices & indices_out, Index & cell_index_out) const {
			Params params;
			if( ! this->value_to_params( v, params, cell_index_out ) ) return false;
			Float scale = Float(ONE<<resl);
			for(size_t i = 0; i < DIM; ++i)	indices_out[i] = static_cast<Index>(params[i]*scale);
			return true;
		}

		void get_indicies_unitcell(Value const & v, Index resl, Indices & indices_out) const {
			Params params;
			this->value_to_params_unitcell( v, params );
			Float scale = Float(ONE<<resl);
			for(size_t i = 0; i < DIM; ++i){
				// this crazy add/subtract avoids round towards 0 so params < 0 behave correctly
				indices_out[i] = static_cast<Index>(params[i]*scale+65536) - 65536;
			}
		}

		Index get_index(Indices const & indices, Index cell_index, Index resl) const {
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i)	index |= util::dilate<DIM>(indices[i]) << i;
			index = index | (cell_index << (DIM*resl));
			return index;
		}

		Index get_index(Value const & v, Index resl) const {
			Index cell_index;
			Indices indices;
			if( !get_indicies(v,resl,indices,cell_index) ) return std::numeric_limits<Index>::max();
			return get_index(indices,cell_index,resl);
		}

		template<class OutIter>
		void push_index(SignedIndices const & indices, Index cell_index, Index resl, OutIter out) const {
			Index index = 0;
			for(size_t i = 0; i < DIM; ++i)	index |= util::dilate<DIM>(indices[i]) << i;
			index = index | (cell_index << (DIM*resl));
			*out = index;
			++out;
		}

		template<class OutIter>
		void get_neighbors(Indices const & indices, Index cell_index, Index resl, OutIter out)  {
			// std::cout << indices.transpose() << std::endl;
			SignedIndices lb = ((indices.template cast<int>()-1).max(     0     ));
			SignedIndices ub = ((indices.template cast<int>()+1).min((1<<resl)-1));
			// std::cout << "LB " << lb.transpose() << std::endl;
			// std::cout << "UB " << ub.transpose() << std::endl;			
			boost::function<void(SignedIndices)> functor;
			functor = boost::bind( & ThisType::template push_index<OutIter>, this, _1, cell_index, resl, out );
			util::NESTED_FOR<DIM>(lb,ub,functor);
		};

		template<class OutIter>
		bool get_neighbors(Value const & v, Index resl, OutIter out)  {
			Indices indices;
			Index cell_index;
			if( !get_indicies(v,resl,indices,cell_index) ) return false;
			get_neighbors(indices,cell_index,resl,out);
			return true;
		};

		template<class OutIter>
		void get_neighbors_unitcell(Value const & v, Index resl, OutIter out)  {
			Indices indices;
			get_indicies_unitcell(v,resl,indices);
			get_neighbors(indices,0,resl,out);
		};

		///////////////////////////////////////
		//// virtual interface functions
		///////////////////////////////////////

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool virtual_set_state(Index index, Index resl) { return set_state(index,resl); }
		///@brief virtual virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@return false iff invalid index
		virtual bool virtual_set_state(
			std::vector<size_t> const & hindices,
			Index cell_index,
			size_t & iindex,
			Index resl
		) {
			Float scale = 1.0 / Float(ONE<<resl);
			Params params;
			for(size_t i = 0; i < DIM; ++i) params[i] = (static_cast<Float>(hindices[iindex]) + 0.5 ) * scale;
			iindex += DIM;
			return this->params_to_value( params, cell_index, this->nonconst_value() );
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
		class Index,
		class Float
	>
	struct NEST<0,Value,ParamMap,StoragePolicy,Index,Float> : 
	              public NestBase<Index>,
	              public ParamMap<0,Value,Index,Float>, 
	              public StoragePolicy<Value>
	{
		typedef Eigen::Array<Float,0,1> Params;		
		NEST(Index cell_size=1) : NestBase<Index>(cell_size) {}
		NEST(std::vector<Value> const & _choices) : 
			NestBase<Index>(_choices.size()),
			ParamMap<0,Value,Index,Float>(_choices)
		{}
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
		//// virtual interface functions
		///////////////////////////////////////

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool virtual_set_state(Index index, Index resl) { return set_state(index,resl); }

		///@brief virtual virtual function to set the state of this nest
		///@detail will consume no indices from hindices vector, using only cell_index
		///        for use in composite data structures containing NestBases
		virtual
		bool
		virtual_set_state(
			std::vector<size_t> const & /*hindices*/,
			Index cell_index,
			size_t & /*iindex*/,
			Index /*resl*/
		) {
			Params params;
			return this->params_to_value( params, cell_index, this->nonconst_value() );
		}
		///@brief return size(resl)
		///@return cell_size for these type
		virtual Index virtual_size(Index /*resl*/) const { return this->cell_size(); }
		///@brief get dimension of this nest
		///@return always 0 for these types
		virtual size_t virtual_dim() const { return 0; }
	};


}
}

#endif
