// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_util_hh
#define INCLUDED_riflib_rifdock_subroutines_util_hh


#include <riflib/types.hh>

#include <utility/io/ozstream.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <numeric/xyzTransform.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <scheme/kinematics/Director.hh>
#include <scheme/search/HackPack.hh>
#include <riflib/util.hh>
#include <riflib/rifdock_typedefs.hh>

#include <rif_dock_test.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/RifFactory.hh>



using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

namespace devel {
namespace scheme {











}}



#endif