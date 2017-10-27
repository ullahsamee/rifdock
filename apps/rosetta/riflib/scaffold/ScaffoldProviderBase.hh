// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_ScaffoldProviderBase_hh
#define INCLUDED_riflib_scaffold_ScaffoldProviderBase_hh

#include <riflib/types.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>



namespace devel {
namespace scheme {

template<
    class _Index
>
struct ScaffoldProviderBase {
    typedef _Index Index;
    ScaffoldProviderBase() {}


    virtual void get_scaffold(Index i) = 0;


};

template<typename _Index>
using ScaffoldProviderOP = shared_ptr<ScaffoldProviderBase< _Index > >;
template<typename _Index>
using ScaffoldProviderCOP = shared_ptr<ScaffoldProviderBase< _Index > const >;


// Key Assumptions of this class:
// - All unvisited TreeIndexes are equal


struct TreeIndex {
    uint64_t depth;
    uint64_t member;
};

struct TreeRelation {
    uint64_t depth;
    uint64_t parent_member;
    uint64_t first_child;
    uint64_t last_child;
}

struct TreeScaffoldProvider :
    public ScaffoldProviderBase<TreeIndex> {

    TreeScaffoldProvider() {
        clear_new_children();
    }


    virtual void get_scaffold(TreeIndex i) = 0;

    virtual void fill_children(TreeIndex i) = 0;

    shared_ptr<std::vector<TreeIndex const> > get_new_children() {
        shared_ptr<std::vector<TreeIndex const> > to_ret = new_children_;
        clear_new_children();
        return to_ret;
        
    }

    shared_ptr<std::vector<TreeIndex const> > new_children_;

protected:
    void add_new_child(TreeIndex const index) {
        new_children_.push_back(index);
    }

private:
    void clear_new_children() {
        new_children_ = shared_ptr<std::vector<TreeIndex const> >( new std::vector<TreeIndex const>() );
    }

};


template<typename _Index>
using ScaffoldProviderOP = shared_ptr<ScaffoldProviderBase< _Index > >;
template<typename _Index>
using ScaffoldProviderCOP = shared_ptr<ScaffoldProviderBase< _Index > const >;

}}

#endif
