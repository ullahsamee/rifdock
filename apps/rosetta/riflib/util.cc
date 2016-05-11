// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#include <riflib/util.hh>
#include <ObjexxFCL/format.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/sasa.hh>
#include <utility/file/file_sys_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

namespace devel {
namespace scheme {

void
add_res_if_not_CGP(
	core::pose::Pose const & pose,
	int ir,
	utility::vector1<core::Size> & res,
	bool verbose = false
){
	if( pose.residue(ir).name3() == "CYD" ){ if(verbose) std::cout << "get_res skip CYD" << std::endl; return; }
	if( pose.residue(ir).name3() == "GLY" ){ if(verbose) std::cout << "get_res skip GLY" << std::endl; return; }
	if( pose.residue(ir).name3() == "PRO" ){ if(verbose) std::cout << "get_res skip PRO" << std::endl; return; }
	res.push_back( ir );
}

utility::vector1<core::Size>
get_res(
	std::string fname,
	core::pose::Pose const & pose,
	bool nocgp
 ){
 	using core::Size;
	utility::vector1<core::Size> res;
	if(fname==""){
		for(Size ir = 1; ir <= pose.n_residue(); ++ir){
			if( nocgp )	add_res_if_not_CGP( pose, ir, res );
			else res.push_back(ir);
		}
		return res;
	}
	runtime_assert_msg( utility::file::file_exists( fname ), "res list file doesn't exist: " + fname );
	utility::io::izstream in(fname);
	std::string s;
	while(in>>s){
		try {
			utility::vector1<std::string> splt = utility::string_split(s,'-');
			// cout << splt << endl;
			if(splt.size()==1){
				// std::cout << "'" << splt.front() << "'" << std::endl;
				Size ir = boost::lexical_cast<Size>(splt.front());
				runtime_assert_msg( ir > 0 && ir <= pose.n_residue(), "residue number out of bounds "+boost::lexical_cast<std::string>(ir) );
				if( nocgp )	add_res_if_not_CGP( pose, ir, res );
				else res.push_back(ir);
			}
			else if(splt.size()==2){
				Size lb = boost::lexical_cast<Size>(splt.front());
				Size ub = boost::lexical_cast<Size>(splt.back());
				if( ub < lb ) std::swap(lb,ub);
				for(Size ir = lb; ir <= ub; ++ir){
					runtime_assert_msg( ir > 0 && ir <= pose.n_residue(), "residue number out of bounds "+boost::lexical_cast<std::string>(ir) );
				if( nocgp )	add_res_if_not_CGP( pose, ir, res );
				else res.push_back(ir);
				}
			} else {
				utility_exit_with_message("get_res: can't parse res line "+s);
			}
		} catch( boost::bad_lexical_cast ){
			utility_exit_with_message("get_res: can't parse res line "+s);
		}
	}
	in.close();

	std::sort( res.begin(), res.end() );
	for( int iri = 1; iri <= res.size(); ++iri ){
		int ir = res[iri];
		// std::cout << "selected res " << ir << " " << pose.residue(ir).name3() << std::endl;
		runtime_assert( pose.residue(ir).name3() != "CYD" );
		runtime_assert( pose.residue(ir).name3() != "GLY" );
		runtime_assert( pose.residue(ir).name3() != "PRO" );
	}
	// utility_exit_with_message("oriestnd");

	return res;
 }

utility::vector1<core::Size> get_res_by_sasa(
	  core::pose::Pose pose
	, bool noloops
	, bool nocgp
 ){
 	core::scoring::dssp::Dssp dssp( pose );
 	dssp.insert_ss_into_pose( pose );
	utility::vector1<core::Size> res;
	// Real calc_per_atom_sasa( pose::Pose const & pose, id::AtomID_Map< Real > & atom_sasa, utility::vector1< Real > & rsd_sasa,
	// Real const probe_radius, bool const use_big_polar_H = false );
	core::id::AtomID_Map< core::Real > atom_sasa;
	utility::vector1< core::Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 2.0 );
	for( int ir = 1; ir <= pose.n_residue(); ++ir ){
		// std::cout << pose.secstruct(ir) << std::endl;
		if( noloops && pose.secstruct(ir) == 'L' ) continue;
		core::Real scsasa = 0;
		for( int ia = 5; ia <= pose.residue(ir).natoms(); ++ia ){
			scsasa += atom_sasa[core::id::AtomID(ia,ir)];
		}
		if( scsasa > 0.0 ){
			if( nocgp )	add_res_if_not_CGP( pose, ir, res );
			else res.push_back(ir);
		} else {
			// std::cout << "pruned res by sasa " << ir << std::endl;
		}
	}
	runtime_assert_msg( res.size() > 0, "no residues selected!" );
	std::cout << "get_res_by_sasa keep " << res.size() << " of " << pose.n_residue() << std::endl;
	std::cout << "select res_by_sasa=resi " << res[1];
	for(int i = 2; i <= res.size(); ++i ) std::cout << "+"<<res[i];
	std::cout << std::endl;

	// utility_exit_with_message("foo");

	return res;
}


std::string KMGT(double const & x, int const & w, int const & d){
	using ObjexxFCL::format::F;
	if( x < 1e3  ) return F( w, d, x/1e0  )+" ";
	if( x < 1e6  ) return F( w, d, x/1e3  )+"K";
	if( x < 1e9  ) return F( w, d, x/1e6  )+"M";
	if( x < 1e12 ) return F( w, d, x/1e9  )+"G";
	if( x < 1e15 ) return F( w, d, x/1e12 )+"T";
	if( x < 1e18 ) return F( w, d, x/1e15 )+"P";
	if( x < 1e21 ) return F( w, d, x/1e18 )+"E";
	if( x < 1e24 ) return F( w, d, x/1e21 )+"Z";
	else           return F( w, d, x/1e24 )+"Y";

}

void pose_to_ala( core::pose::Pose & pose ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA") );
	for( int ir = 1; ir <= pose.n_residue(); ++ir ){
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *ala, true );
	}
}
void pose_to_ala( core::pose::Pose & pose, utility::vector1<core::Size> const & res_sel ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA") );
	for( int iri = 1; iri <= res_sel.size(); ++iri ){
		int ir = res_sel[iri];
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *ala, true );
	}
}


std::string
open_for_read_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::izstream & in
){
	for( auto const & dir : path ){
		if( utility::file::file_exists( dir+"/"+fname ) ){
			in.close();
			in.clear();
			in.open( dir+"/"+fname, std::ios::binary );
			if( in.good() ){
				// std::cout << "open_for_read_on_path: opened " << dir + "/" + fname << std::endl;
				return dir+"/"+fname;
			}
		}
	}
	return std::string();
}

std::string
open_for_write_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::ozstream & out,
	bool create_directorys /* = false */
){
	for( auto const & dir: path ){
		if( create_directorys && !utility::file::file_exists( dir ) ){
			utility::file::create_directory_recursive( dir );
		}
		if( utility::file::file_exists( dir ) ){
			out.close();
			out.clear();
			out.open( dir+"/"+fname, std::ios::binary );
			if( out.good() ){
				// std::cout << "open_for_write_on_path: opened " << dir + "/" + fname << std::endl;
				return dir+"/"+fname;
			}
		}
	}
	return std::string();

}









}
}


