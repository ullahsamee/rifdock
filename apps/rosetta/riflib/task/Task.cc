// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/task/Task.hh>

#include <scheme/util/type.hh>

#include <utility/string_util.hh>


namespace devel {
namespace scheme {


std::string
Task::name() const {
    std::string long_name = ::scheme::type(*this);

    utility::vector1<std::string> split = utility::string_split(long_name, ':');
    std::string short_name = split[split.size()];
    return short_name;
}




}}
