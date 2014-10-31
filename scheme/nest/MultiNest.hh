#ifndef INCLUDED_scheme_nest_MultiNest_HH
#define INCLUDED_scheme_nest_MultiNest_HH

#include "scheme/nest/NEST.hh"

#include <boost/assert.hpp>

namespace scheme { namespace nest {



	///@brief Base class for NEST
	///@tparam Index index type
	///@detail Base class for NEST, virtual NEST interface
	template<
		class Index=uint64_t,
		class BigIndex=uint64_t
	>
	struct MultiNest : public NestBase<Index> {

		typedef NestBase<Index> Nest;
		typedef shared_ptr<Nest> Nestp;
		typedef std::vector<Nestp> Nests;
		typedef std::vector<Index> Indices;

		Index dim_, tot_ncells_;
		Nests nests_;
		Indices ncells_, ncells_pref_sum_;

		MultiNest(){ init(); }
		MultiNest(Nests const & nests){ init(nests); }

		void add_nest( Nestp nestp ){
			nests_.push_back(nestp);
			init();
		}

		void init( Nests const & nests ){ nests_ = nests; init(); }

		void init(){
			ncells_.resize(nests_.size());
			ncells_pref_sum_.resize(nests_.size());
			dim_ = 0;
            if(nests_.size()==0) return;
			tot_ncells_ = 1;
			for(size_t i = 0; i < nests_.size(); ++i){
				Nestp nest = nests_[i];
				dim_ += nest->virtual_dim();
				ncells_[i] = nest->virtual_num_cells();
				tot_ncells_ *= ncells_[i];
			}
			ncells_pref_sum_[0] = 1;
			for(size_t i = 1; i < nests_.size(); ++i){
				ncells_pref_sum_[i] = ncells_pref_sum_[i-1] * ncells_[i-1];
			}
		}

		template<class Anys>
		void get_states(BigIndex const & index, Index resl, Anys & anys) const {
			BOOST_VERIFY(anys.size()==nests_.size());
			BOOST_VERIFY( index < size(resl) );
			BigIndex big_cell_index = index >> dim_*resl;
			BigIndex const one(1);
			BigIndex hier_index = index & ((one<<(dim_*resl))-one);
			// std::cout << "cell_index: " << big_cell_index << std::endl;
			// std::cout << "hier_index: " << hier_index << std::endl;	
			Indices hier_indices;
			expand_index( hier_index, hier_indices );
			size_t iindex = 0;
			for(size_t i = 0; i < nests_.size(); ++i){
				Index cell_index = big_cell_index / ncells_pref_sum_[i] % ncells_[i];
				// std::cout << "nest " << i << " " << cell_index << std::endl;
				nests_[i]->virtual_get_state( hier_indices, cell_index, iindex, resl, anys[i] );
			}
		}

		void expand_index( BigIndex const & index, Indices & out ) const {
			out.resize(dim_);
			for(Index i = 0; i < dim_; ++i){
				out[i] = util::undilate( dim_, index>>i );
			}
		 }

		///@brief need virtual destructor
		virtual ~MultiNest(){}

		///@brief virtual virtual function to set the state of this nest
		///@returns false if invalid index
		virtual bool 
		virtual_get_state(
			Index index, 
			Index resl, 
			boost::any & result
		){
			std::cout << "not implemented yet" << std::endl;
			std::exit(-1);
		}

		///@brief virtual virtual function to set the state of this nest
		///@detail will consume DIM indices from hindices vector, starting at iindex, then will increment iindex
		///        for use in composite data structures containing NestBases
		///@returns false if invalid index
		virtual	bool 
		virtual_get_state(
			std::vector<Index> const & indices,
			Index cell_index,
			size_t & iindex,
			Index resl,
			boost::any & result
		){
			std::cout << "not implemented yet" << std::endl;
			std::exit(-1);
		}

		///@brief get the total size of this NEST at resl
		///@return number of possible states at depth resl
		virtual Index 
		virtual_size(Index resl) const { 
			Index s = (Index)size(resl);
			BOOST_VERIFY( s <= BigIndex( std::numeric_limits<Index>::max() ) );
			return s;
		}

		///@brief get the number of cells in this nest
		virtual Index 
		virtual_num_cells() const {
			return tot_ncells_;
		}

		///@brief get the dimension of this nest
		///@return dimension of Nest
		virtual size_t 
		virtual_dim() const { return dim(); }
		
		size_t dim() const { return dim_; }
		BigIndex size( Index resl ) const {
			return BigIndex((unsigned long)tot_ncells_) * ((BigIndex(1))<<(resl*dim_));
		}
	};


// uint64_t const BOGUS_INDEX = 0xFFFFFFFFFFFFFFFFul;

// using std::cout;
// using std::endl;
// using ObjexxFCL::format::F;

// static const uint64_t ONE(1);
// static const uint64_t TWO(2);
// static const uint64_t THREE(3);




// void CompositeNest::set_position( NestIndex const & index, uint64_t resl ) {
// 	uint64_t bindex0 = (index >> (dim_*resl)).touint64();
// 	NestIndex hindex0 = index & (NestIndex(1<< (dim_*resl))-1);
// 	vector<uint64_t> hindex;
// 	expand_index(hindex0,dim_,hindex);
// 	// cout << index << " B" << bindex0 << " H" << hindex0;
// 	// BOOST_FOREACH(uint64_t i,hindex) cout << " " << i;
// 	// cout << endl;
// 	uint64_t nb0 = 1;
// 	vector<uint64_t>::iterator ih = hindex.begin();
// 	vector<uint64_t> atomicindex;
// 	for(unsigned ibase = 0; ibase < nbase_.size(); ++ibase){
// 		uint64_t const & nb = nbase_[ibase];
// 		uint64_t const thisbindex = (( bindex0 / nb0 ) % nb);
// 		// cout << "basei " << (( bindex0 / nb0 ) % nb) << endl;
// 		if(is_dependent_[ibase]){
// 			; // pass
// 		} else if(is0dim_[ibase]){
// 			atomicindex.push_back(thisbindex);
// 		} else {
// 			atomicindex.push_back( *ih + (thisbindex << resl) );
// 			++ih;
// 		}
// 		nb0 *= nb;
// 	}
// 	// cout << " resulting indices from " << index << ": ";
// 	// BOOST_FOREACH(uint64_t i,atomicindex) cout << " " << i;
// 	// cout << endl;
// 	int tmp=0;
// 	BOOST_FOREACH(NestOP gop,nests_) gop->set_position(atomicindex,tmp,resl);

// 	runtime_assert(tmp==(int)atomicindex.size());
// 	BOOST_FOREACH(Size2 ds,dependent_map_){
// 		DependentNest * dNest = dynamic_cast<DependentNest*>(nests_[ds.first ]());
// 		runtime_assert( dNest != 0 );
// 		dNest->modify_xform();
// 	}
// 	// vector<NestOP>::const_iterator ig = nests_.begin();
// 	// BOOST_FOREACH(Xform x,out) cout << x.t << " " << (*(ig++))->name() <<  endl;
//  }



}}

#endif
