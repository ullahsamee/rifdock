// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_RifGeneratorSimpleHbonds_hh
#define INCLUDED_riflib_rif_RifGeneratorSimpleHbonds_hh


#include <riflib/rif/RifGenerator.hh>
#include <riflib/rif/make_hbond_geometries.hh>

namespace devel {
namespace scheme {
namespace rif {

struct RifGeneratorSimpleHbondsOpts {
	float tip_tol_deg = 30.0;
	float rot_samp_resl = 6.0;
	float rot_samp_range = 360.0;
	float hbond_cart_sample_hack_range = 0.25;
	float hbond_cart_sample_hack_resl = 0.25;
	float score_threshold = 0.0;
	float dump_fraction = 0.0;
	float hbond_weight = 2.0;
	float upweight_multi_hbond = 0;
	float min_hb_quality_for_satisfaction = -0.6;
	float long_hbond_fudge_distance = 0.0;
	bool debug = false;
	bool dump_bindentate_hbonds = false;
	bool report_aa_count = false;
	int hbgeom_max_cache = -1;
};

struct HBJob {
	std::string don, acc, don_or_acc;
	int nrots;
	int ires;
    std::string hbgeomtag;
	bool operator<( HBJob const & other ) const { return nrots > other.nrots; }
};


struct RifGeneratorSimpleHbonds : public RifGenerator {

	utility::vector1<std::string> donresn_user;
	utility::vector1<std::string> accresn_user;
	RifGeneratorSimpleHbondsOpts opts;


	RifGeneratorSimpleHbonds(
		  utility::vector1<std::string> donresn_user
		, utility::vector1<std::string> accresn_user
		, RifGeneratorSimpleHbondsOpts opts
	)
	: donresn_user( donresn_user )
	, accresn_user( accresn_user )

	, opts( opts )
	{}

	void
	prepare_hbgeoms( 
	    std::vector<HBJob> const & hb_jobs,
	    int start_job,
	    int end_job,    // python style numbering. To do all jobs, specify 0, hb_jobs.size()
	    std::map< std::string, utility::vector1< RelRotPos > * > & hbond_geoms_cache,
	    std::map< std::string, omp_lock_t > & hbond_io_locks,
	    omp_lock_t & cout_lock,
	    omp_lock_t & io_lock,
	    omp_lock_t & pose_lock,
	    omp_lock_t & hacky_rpms_lock,
	    omp_lock_t & hbond_geoms_cache_lock,
	    RifGenParamsP const & params
	);

	void generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	) override;

};

}
}
}

#endif
