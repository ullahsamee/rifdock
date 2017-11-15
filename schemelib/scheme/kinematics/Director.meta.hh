#ifndef INCLUDED_kinematics_Director_meta.HH
#define INCLUDED_kinematics_Director_meta.HH

#include <scheme/util/meta/util.hh>
template <class __Director>
using _DirectorPosition = typename ::scheme::util::meta::remove_pointer<__Director>::type::Position;
template <class __Director>
using _DirectorIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::Index;
template <class __Director>
using _DirectorBigIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::BigIndex;
template <class __Director>
using _DirectorScaffoldIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::ScaffoldIndex;
#endif
