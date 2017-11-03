// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_scaffold_ScaffoldProviderBase_hh
#define INCLUDED_scaffold_ScaffoldProviderBase_hh

#include <scheme/types.hh>

#include <string>
#include <vector>
#include <limits>
#include <boost/any.hpp>



namespace scheme {
namespace scaffold {

const uint64_t BOGUS_INDEX = std::numeric_limits<uint64_t>::max();

template<
    class _Scaffold,
    class _ScaffoldIndex,
    class _ScaffoldIndexLimits
>
struct ScaffoldProviderBase {
    typedef _Scaffold Scaffold;
    typedef _ScaffoldIndex ScaffoldIndex;
    typedef _ScaffoldIndexLimits ScaffoldIndexLimits;
    ScaffoldProviderBase() {}


    virtual shared_ptr<Scaffold const> get_scaffold(ScaffoldIndex i) = 0;

    virtual ScaffoldIndexLimits get_scaffold_index_limits() = 0;

};

template<typename _Scaffold, typename _ScaffoldIndex, typename _ScaffoldIndexLimits>
using ScaffoldProviderOP = shared_ptr<ScaffoldProviderBase< _Scaffold, _ScaffoldIndex, _ScaffoldIndexLimits > >;
template<typename _Scaffold, typename _ScaffoldIndex, typename _ScaffoldIndexLimits>
using ScaffoldProviderCOP = shared_ptr<ScaffoldProviderBase< _Scaffold, _ScaffoldIndex, _ScaffoldIndexLimits > const >;


// Key Assumptions of this class:


struct TreeIndex {
    uint64_t depth;
    uint64_t member;
    TreeIndex() : depth(BOGUS_INDEX), member(BOGUS_INDEX) {};
};

struct TreeRelation {
    uint64_t depth;
    uint64_t parent_member;
    uint64_t first_child;
    uint64_t last_child;
};

typedef std::vector<Bounds<uint64_t>> TreeLimits;


template<
    class _Scaffold
>
struct TreeScaffoldProvider :
    public ScaffoldProviderBase<_Scaffold, TreeIndex, TreeLimits> {
    typedef _Scaffold Scaffold;


    virtual shared_ptr<Scaffold const> get_scaffold(TreeIndex i) = 0;

    virtual void fill_children(TreeIndex i) = 0;

    virtual TreeLimits get_scaffold_index_limits() = 0;


};


template<typename _Scaffold>
using TreeProviderOP = shared_ptr<TreeScaffoldProvider< _Scaffold> >;
template<typename _Scaffold>
using TreeProviderCOP = shared_ptr<TreeScaffoldProvider< _Scaffold> const >;

inline
uint64_t scaffold_index_default_value(uint64_t) {
    return BOGUS_INDEX;
}

inline
TreeIndex scaffold_index_default_value(TreeIndex) {
    return TreeIndex();
}




}}

#endif
