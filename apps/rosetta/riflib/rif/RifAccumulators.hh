// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_RifAccumulators_hh
#define INCLUDED_riflib_rif_RifAccumulators_hh

#include <riflib/rif/RifGenerator.hh>
#include <riflib/RifFactory.hh>

namespace devel {
namespace scheme {
namespace rif {


template<class XMap>
struct RIFAccumulatorMapThreaded : public RifAccumulator {

	typedef typename XMap::Map Map;
	shared_ptr<RifFactory const> rif_factory_;
	std::vector< Map > to_insert_;
	std::vector<int64_t> nsamp_;
	float scratch_size_M_;
	uint64_t N_motifs_found_;

	shared_ptr<XMap> xmap_ptr_;

	RIFAccumulatorMapThreaded(
		shared_ptr<RifFactory const> rif_factory,
		float cart_resl,
		float ang_resl,
		float cart_bound,
		size_t scratch_size_M=8000
	)
		: rif_factory_(rif_factory)
		, scratch_size_M_(scratch_size_M)
	 	, N_motifs_found_(0)
	{
		clear();
		xmap_ptr_ = make_shared<XMap>( cart_resl, ang_resl );
	}

	uint64_t n_motifs_found() const override { return N_motifs_found_ + total_samples(); }

	shared_ptr<RifBase> rif() const override {
		shared_ptr<RifBase> r = rif_factory_->create_rif();
		r->set_xmap_ptr( xmap_ptr_ );
		return r;
	}

	void insert( devel::scheme::EigenXform const & x, float score, int32_t rot, int sat1, int sat2, bool single_thread ) override {
		if( score > 0.0 ) return;
		uint64_t const key = xmap_ptr_->hasher_.get_key( x );
		typename XMap::Map & map_for_this_thread( single_thread ? xmap_ptr_->map_ : to_insert_[ omp_get_thread_num() ] );
		// std::cerr << "INSERT mapsize: " << map_for_this_thread.size() << " thread: " << omp_get_thread_num() << " nmaps: " << to_insert_.size() << std::endl;
		typename XMap::Map::iterator iter = map_for_this_thread.find(key);
		if( iter == map_for_this_thread.end() ){
			typename XMap::Value value;
			value.add_rotamer( rot, score, sat1, sat2 );


         //    std::vector<int> sats;
        	// value.rotamer_sat_groups( rot, sats );
        	// std::stringstream ss;
        	// ss << "sats: "
        	// for (auto val : sats) {
        	// 	ss << " : " << val;
        	// }
        	// ss << std::endl;
        	// std::cout << ss.str();
			// runtime_assert( value.score(0) > -9.0 );
			map_for_this_thread.insert( std::make_pair( key, value ) );
		} else {
			// for( int i = 0; i < iter->second.maxsize(); ++i ) runtime_assert( iter->second.score(i) > -9.0 );
			// std::cout << score << std::endl;
			iter->second.add_rotamer( rot, score, sat1, sat2 );
			// for( int i = 0; i < iter->second.maxsize(); ++i ) runtime_assert( iter->second.score(i) > -9.0 );
		}
		++nsamp_[ omp_get_thread_num() ];
	}

	int64_t total_samples() const override {
		int64_t tot = 0;
		for( int i = 0; i < nsamp_.size(); ++i ) tot += nsamp_[i];
		return tot;
	}

	bool need_to_condense() const override {
		return mem_use() > uint64_t(scratch_size_M_)*uint64_t(1024*1024);
	}

	void condense() override {
		using ObjexxFCL::format::I;
		for( int i = 0; i < to_insert_.size(); ++i ){
			// std::cout << I(3,i+1) << " of " << to_insert_.size() << " progress: ";
			int64_t const out_interval = std::max<int64_t>(1,to_insert_[i].size()/100);
			int64_t count = 0;
			BOOST_FOREACH( typename XMap::Map::value_type const & value, to_insert_[i] ){
				// if( ++count%out_interval == 0 ){ std::cout << '*'; std::cout.flush(); }
				typename XMap::Key const key = value.first;
				typename XMap::Value const & rotsc = value.second;
				typename XMap::Map::iterator iter = xmap_ptr_->map_.find(key);
				if( iter == xmap_ptr_->map_.end() ){
					xmap_ptr_->map_.insert( std::make_pair( key, rotsc ) );
					// for( int i = 0; i < rotsc.maxsize(); ++i ) runtime_assert( rotsc.score(i) > -9.0 );
				} else {
					iter->second.merge( rotsc );
					// for( int i = 0; i < iter->second.maxsize(); ++i ) runtime_assert( iter->second.score(i) > -9.0 );
				}

			}
			// std::cout << std::endl;
		}
	}

	void report( std::ostream & out ) const override {
		out << "RIFAccum nrots: " << devel::scheme::KMGT(n_motifs_found())
		    << " mem: " << devel::scheme::KMGT(mem_use())
		    << " rif_mem: " << devel::scheme::KMGT(xmap_ptr_->mem_use()) << std::endl;
	}

	// inclusive on the ranges
	uint64_t count_these_irots( int irot_low, int irot_high ) const {
		uint64_t count = 0;
		for ( auto pair : xmap_ptr_->map_ ) {
			count += pair.second.count_these_irots( irot_low, irot_high );
		}
		return count;
	}

	// This isn't threadsafe at all!!!
	std::set<size_t> get_sats_of_this_irot( devel::scheme::EigenXform const & x, int irot ) const override {

		std::set<size_t> sats;

		uint64_t const key = xmap_ptr_->hasher_.get_key( x );

		typename XMap::Map & map_for_this_thread(xmap_ptr_->map_);
		typename XMap::Map::iterator iter = map_for_this_thread.find(key);

		typedef typename XMap::Value Value;


		if( iter == map_for_this_thread.end() ){

			return sats;
		} else {

			Value const & rotscores = iter->second;
			static int const Nrots = Value::N;
			// typedef typename Value::RotScore RotScore;

			for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
				if( rotscores.empty(i_rs) ) {
					break;
				}

                int _irot = rotscores.rotamer(i_rs);
                if (_irot != irot) continue;

				sats.insert(255);

				std::vector<int> sat_groups;
				rotscores.rotamer_sat_groups( i_rs, sat_groups );

				for ( int number : sat_groups ) {
					sats.insert(number);
				}
			}
		}

		return sats;

	}

	uint64_t mem_use() const {
		uint64_t mem = 0;
		for( int i = 0; i < to_insert_.size(); ++i ){
			mem += to_insert_[i].bucket_count()*(sizeof(typename XMap::Map::value_type));
		}
		return mem;
	}


	void clear() override {
		for( int i = 0; i < to_insert_.size(); ++i ) to_insert_[i].clear();
		to_insert_.clear();
		nsamp_.clear();

		to_insert_.reserve( devel::scheme::omp_max_threads_1() );
		for( int i = 0; i < devel::scheme::omp_max_threads_1(); ++i ){
			Map m;
			m.set_empty_key( std::numeric_limits<uint64_t>::max() );
			to_insert_.push_back(m);
		}
		nsamp_.resize( devel::scheme::omp_max_threads_1(), 0 );
	}

	void checkpoint( std::ostream & out ) override {
		using devel::scheme::KMGT;
		// uint64_t m = mem_use();
		// if( m > uint64_t(scratch_size_M_)*uint64_t(1024*1024) ){ // time to clean up a bit...
			// out << "mem use " << float(m)/1024.0/1024.0 << "M is above threshold " << scratch_size_M_ << "M, time to condense" << std::endl;
			out << '<'; out.flush();
			condense();
			N_motifs_found_ += total_samples();
			clear();
			out << '>'; out.flush();
			// out << "RIF so far: " << " non0 in RIF: " << KMGT(xmap_ptr_->size()) << " N_motifs_found_: "
			    // << KMGT(N_motifs_found_) << " coverage: " << (double)N_motifs_found_/xmap_ptr_->size()
			    // << ", mem_use: " << KMGT(xmap_ptr_->mem_use())
			    // << ", load: " << KMGT(xmap_ptr_->map_.size()*1.f/xmap_ptr_->map_.bucket_count())
			    // << std::endl;
		// }
		// report( out );
	}

};



}
}
}

#endif
