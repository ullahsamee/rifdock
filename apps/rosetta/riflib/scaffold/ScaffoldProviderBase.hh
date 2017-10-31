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
    class _Scaffold,
    class _ScaffoldCache,
    class _Index,
    class _IndexLimits
>
struct ScaffoldProviderBase {
    typedef _Scaffold Scaffold;
    typedef _ScaffoldCache ScaffoldCache;
    typedef _Index Index;
    typedef _IndexLimits IndexLimits;
    ScaffoldProviderBase() {}


    virtual Scaffold get_scaffold(Index i) = 0;
    virtual ScaffoldCache get_scaffold_cache() = 0;

    virtual IndexLimits get_index_limits() = 0;

};

template<typename _Scaffold, typename _ScaffoldCache, typename _Index, typename _IndexLimits>
using ScaffoldProviderOP = shared_ptr<ScaffoldProviderBase< _Scaffold, _ScaffoldCache, _Index, _IndexLimits > >;
template<typename _Scaffold, typename _ScaffoldCache, typename _Index, typename _IndexLimits>
using ScaffoldProviderCOP = shared_ptr<ScaffoldProviderBase< _Scaffold, _ScaffoldCache, _Index, _IndexLimits > const >;


// Key Assumptions of this class:


struct TreeIndex {
    uint64_t depth;
    uint64_t member;
};

struct TreeRelation {
    uint64_t depth;
    uint64_t parent_member;
    uint64_t first_child;
    uint64_t last_child;
};

typedef std::vector<Range> TreeLimits;


template<
    class _Scaffold,
    class _ScaffoldCache
>
struct TreeScaffoldProvider :
    public ScaffoldProviderBase<_Scaffold, _ScaffoldCache, TreeIndex, TreeLimits> {
    typedef _Scaffold Scaffold;
    typedef _ScaffoldCache ScaffoldCache;


    virtual Scaffold get_scaffold(TreeIndex i) = 0;
    virtual ScaffoldCache get_scaffold_cache() = 0;

    virtual void fill_children(TreeIndex i) = 0;

    virtual TreeLimits get_index_limits() = 0;


};


template<typename _Scaffold, typename _ScaffoldCache>
using TreeProviderOP = shared_ptr<TreeScaffoldProvider< _Scaffold, _ScaffoldCache > >;
template<typename _Scaffold, typename _ScaffoldCache>
using TreeProviderCOP = shared_ptr<TreeScaffoldProvider< _Scaffold, _ScaffoldCache > const >;



struct ScaffoldDataCache {
    std::vector<std::vector<float> > const * rotamer_energies_1b_;
    
};









}}

#endif
