// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:



#include <riflib/seeding_util.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <core/import_pose/import_pose.hh>


namespace devel {
namespace scheme {




shared_ptr<std::vector<EigenXform>>
setup_seeding_positions( RifDockOpt & opt, ProtocolData & pd, ScaffoldProviderOP & scaffold_provider, int iscaff ) {

    shared_ptr<std::vector<EigenXform>> seeding_positions = make_shared<std::vector<EigenXform>>();

    if ( opt.seed_with_these_pdbs.size() > 0 ) {
        if ( opt.seed_include_input ) {
            std::cout << "Adding seeding position: " << "_SP_input" << std::endl;
            seeding_positions->push_back(EigenXform::Identity());
            pd.seeding_tags.push_back("_input");
        }

        core::pose::Pose const & input = *scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scaffold_unmodified_p;

        for ( std::string const & pdb : opt.seed_with_these_pdbs ) {
            std::string tag = "_SP_" + pdb_name( pdb );
            std::cout << "Adding seeding position: " << tag << std::endl;
            core::pose::Pose pose;
            core::import_pose::pose_from_file( pose, pdb );
            std::cout << pose.size() << " " << input.size() << std::endl;
            EigenXform xform = find_xform_from_identical_pose_to_pose( input, pose, 1 );
            seeding_positions->push_back(xform);
            pd.seeding_tags.push_back( tag );
        }

        return seeding_positions;
    }


    if( opt.seeding_fnames.size() ){

        std::string seeding_fname = "";
        if( opt.seeding_fnames.size() == opt.scaffold_fnames.size() ) {
            seeding_fname = opt.seeding_fnames.at(iscaff);
        } else {
            utility_exit_with_message( "-seeding_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_seeding_file(seeding_fname, *seeding_positions, opt.seeding_by_patchdock), "Faild to parse the seeding file!!!");
        
        Eigen::Vector3f scaffold_center = scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scaffold_center;

        EigenXform x(EigenXform::Identity());
        x.translation() = scaffold_center;
        for( auto & t : *seeding_positions ) t = t * x;

        return seeding_positions;
    }



    return nullptr;

}




bool 
parse_exhausitive_searching_file(
    std::string fname, 
    std::vector<std::pair< int64_t, devel::scheme::EigenXform > > & searching_positions, 
    double maximum_ang/* = 999*/
) {
        
    // the seeding
        
    if (fname == "" ) {
        return false;
    }
    runtime_assert_msg(utility::file::file_exists( fname ), "exhausitive searching file does not exits: " + fname );
    
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
    std::string s;
    devel::scheme::EigenXform xform;
    while (std::getline(in, s)) {
        if (s.empty()) continue;
            
        utility::vector1<std::string> splt = utility::string_split_simple(s, ' ');
        if (splt[1].find("#") == 0  || splt.size() == 0) {
            continue;
        }
        
        if ( std::abs( utility::string2float(splt[3]) ) > maximum_ang ) {
            continue;
        }
        
        xform.linear().row(0) = Eigen::Vector3f(utility::string2float(splt[4]), utility::string2float(splt[5]), utility::string2float(splt[6]) );
        xform.linear().row(1) = Eigen::Vector3f(utility::string2float(splt[7]), utility::string2float(splt[8]), utility::string2float(splt[9]) );
        xform.linear().row(2) = Eigen::Vector3f(utility::string2float(splt[10]), utility::string2float(splt[11]), utility::string2float(splt[12]) );
        xform.translation() = Eigen::Vector3f(utility::string2float(splt[13]), utility::string2float(splt[14]), utility::string2float(splt[15]) );
        
        searching_positions.push_back( std::pair< int64_t, devel::scheme::EigenXform >( utility::string2int(splt[2]), xform) );
    }
        
    return true;
}

    


bool 
parse_seeding_file(
    std::string fname, 
    std::vector<devel::scheme::EigenXform> & seeding_positions, 
    bool seeding_by_patchdock
) {
        
    // the seeding
    
    if (fname == "" ) {
        return false;
    }
    runtime_assert_msg(utility::file::file_exists( fname ), "seeding file does not exits: " + fname );
    
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
        std::string s;
        devel::scheme::EigenXform xform;
        if ( seeding_by_patchdock ){

                bool flag = false;
                while ( std::getline(in, s)){
                        if (s.empty()) continue;

                        utility::vector1<std::string> splt = utility::quoted_split( s );
                        
                        if (!flag && splt[1].find("#") == 0  && splt[3].find("score") == 0 ) {
                                flag = true;
                                continue;
                        }

                        if(!flag) continue;

                        float cx = cos(utility::string2float(splt[25]));
                        float cy = cos(utility::string2float(splt[26]));
                        float cz = cos(utility::string2float(splt[27]));
                        float sx = sin(utility::string2float(splt[25]));
                        float sy = sin(utility::string2float(splt[26]));
                        float sz = sin(utility::string2float(splt[27]));
                        float tx = utility::string2float(splt[28]);
                        float ty = utility::string2float(splt[29]);
                        float tz = utility::string2float(splt[30]);
                        
                        xform.linear().row(0) = Eigen::Vector3f( cz * cy, -sy * sx * cz - sz * cx, -sy * cx *cz + sz *sx);
                        xform.linear().row(1) = Eigen::Vector3f( sz * cy, -sy * sx * sz + cx * cz, -sy * cx * sz - sx * cz);
                        xform.linear().row(2) = Eigen::Vector3f( sy, cy * sx, cy * cx);
                        xform.translation() = Eigen::Vector3f( tx, ty, tz );

                        seeding_positions.push_back(xform);
                }
        } else {
                // the format of the seeding file is Rosetta Xforms. the easy case
                while (std::getline(in, s)) {
                        if (s.empty()) continue;
                        
                        utility::vector1<std::string> splt = utility::quoted_split( s );
                        if (splt[1].find("#") == 0  || splt.size() == 0) {
                                continue;
                        }
        
                        xform.linear().row(0) = Eigen::Vector3f(utility::string2float(splt[1]), utility::string2float(splt[2]), utility::string2float(splt[3]) );
                        xform.linear().row(1) = Eigen::Vector3f(utility::string2float(splt[4]), utility::string2float(splt[5]), utility::string2float(splt[6]) );
                        xform.linear().row(2) = Eigen::Vector3f(utility::string2float(splt[7]), utility::string2float(splt[8]), utility::string2float(splt[9]) );
                        xform.translation() = Eigen::Vector3f(utility::string2float(splt[10]), utility::string2float(splt[11]), utility::string2float(splt[12]) );
        
                        seeding_positions.push_back(xform);
                }
    }
        
    return true;
}
    





}}






