// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#ifndef INCLUDED_riflib_scaffold_ExtraScaffoldData_hh
#define INCLUDED_riflib_scaffold_ExtraScaffoldData_hh


#include <scheme/types.hh>

#include <riflib/HSearchConstraints.hh>
#include <riflib/AtomsCloseTogetherManager.hh>

#include <vector>


namespace devel {
namespace scheme {

struct ExtraScaffoldData {
    ExtraScaffoldData () :
        csts(),
        force_scaffold_center( Eigen::Vector3f {std::numeric_limits<double>::quiet_NaN(), 0, 0} ),
        per_rotamer_custom_energies_p( nullptr ),
        allowed_irot_at_ires_p( nullptr )
    {}


    std::vector<CstBaseOP> csts;
    Eigen::Vector3f force_scaffold_center;
    std::shared_ptr< std::vector< std::vector<float> > > per_rotamer_custom_energies_p;
    shared_ptr<std::vector<std::vector<bool>>> allowed_irot_at_ires_p;
    shared_ptr<std::vector<bool>> ala_disallowed_p;
    shared_ptr<std::vector<SimpleAtom>> clash_context_p;
    shared_ptr<std::vector<AtomsCloseTogetherManager>> atoms_close_together_managers_p;


    void add_to_clash_context( std::vector<SimpleAtom> const & clash_context ) {
        if ( ! clash_context_p ) {
            clash_context_p = make_shared<std::vector<SimpleAtom>>( clash_context );
        } else {
            clash_context_p->insert(clash_context_p->end(), clash_context.begin(), clash_context.end() );
        }
    }

    void adjust_clash_context( Eigen::Vector3f const & to_move ) {
        if ( ! clash_context_p ) return;

        for ( SimpleAtom & atom : *clash_context_p ) {
            atom.set_position( atom.position() + to_move );
        }
    }


    void accumulate_per_rotamer_custom_energies( std::shared_ptr< std::vector< std::vector<float> > > const & per_rotamer_custom_energies_p_in ) {
        if ( ! per_rotamer_custom_energies_p_in ) return;

        if ( ! per_rotamer_custom_energies_p ) per_rotamer_custom_energies_p = per_rotamer_custom_energies_p_in;

        runtime_assert( per_rotamer_custom_energies_p_in->size() == per_rotamer_custom_energies_p->size() );

        for ( core::Size irow = 0; irow < per_rotamer_custom_energies_p->size(); irow++ ) {
            runtime_assert( per_rotamer_custom_energies_p_in->at(irow).size() == per_rotamer_custom_energies_p->at(irow).size() );
            for ( core::Size icol = 0; icol < per_rotamer_custom_energies_p->at(irow).size(); icol++ ) {
                per_rotamer_custom_energies_p->at(irow)[icol] += per_rotamer_custom_energies_p_in->at(irow)[icol];
            }
        }
    }


    void accumulate_allowed_irot_at_ires( std::shared_ptr< std::vector< std::vector<bool> > > const & allowed_irot_at_ires_p_in ) {
        if ( ! allowed_irot_at_ires_p_in ) return;

        if ( ! allowed_irot_at_ires_p ) allowed_irot_at_ires_p = allowed_irot_at_ires_p_in;

        runtime_assert( allowed_irot_at_ires_p_in->size() == allowed_irot_at_ires_p->size() );

        for ( core::Size irow = 0; irow < allowed_irot_at_ires_p->size(); irow++ ) {
            runtime_assert( allowed_irot_at_ires_p_in->at(irow).size() == allowed_irot_at_ires_p->at(irow).size() );
            for ( core::Size icol = 0; icol < allowed_irot_at_ires_p->at(irow).size(); icol++ ) {
                allowed_irot_at_ires_p->at(irow)[icol] = allowed_irot_at_ires_p->at(irow)[icol] && allowed_irot_at_ires_p_in->at(irow)[icol];
            }
        }
    }

    void accumulate_ala_disallowed( shared_ptr<std::vector<bool>> const & ala_disallowed_p_in ) {
        if ( ! ala_disallowed_p_in ) return;

        if ( ! ala_disallowed_p ) ala_disallowed_p = ala_disallowed_p_in;

        runtime_assert( ala_disallowed_p_in->size() == ala_disallowed_p->size() );

        for ( core::Size irow = 0; irow < ala_disallowed_p->size(); irow++ ) {
            ala_disallowed_p->at(irow) = ala_disallowed_p->at(irow) || ala_disallowed_p_in->at(irow);
        }
    }


};


}
}

#endif