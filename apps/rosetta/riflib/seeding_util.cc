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
            if ( opt.write_seed_to_output ) {
                pd.seeding_tags.push_back("_input");
            }
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
            if ( opt.write_seed_to_output ) {
                pd.seeding_tags.push_back( tag );
            }
        }

        if ( seeding_positions->size() >= 4294967296 ) {
            utility_exit_with_message("Too many seeding positions!!!!");
        }

        return seeding_positions;
    }


    if( opt.seeding_fnames.size() ){

        std::string seeding_fname = "";
        if ( opt.seeding_fnames.size() == 1 ) {
            seeding_fname = opt.seeding_fnames.front();
        } else if( opt.seeding_fnames.size() == opt.scaffold_fnames.size() ) {
            seeding_fname = opt.seeding_fnames.at(iscaff);
        } else {
            utility_exit_with_message( "-seeding_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_seeding_file(seeding_fname, *seeding_positions, opt.seeding_by_patchdock, opt.patchdock_min_sasa, opt.patchdock_top_ranks), "Faild to parse the seeding file!!!");
        
        Eigen::Vector3f scaffold_center = scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scaffold_center;

        if ( ! opt.apply_seeding_xform_after_centering ) {
            EigenXform x(EigenXform::Identity());
            x.translation() = scaffold_center;
            for( auto & t : *seeding_positions ) t = t * x;
        }

        if ( opt.write_seed_to_output ) {
            size_t digits = boost::str(boost::format("%i")%(seeding_positions->size()-1)).length();
            std::string format = "%0" + boost::str(boost::format("%i")%digits) + "i";
            for ( size_t i = 0; i < seeding_positions->size(); i++ ) {
                pd.seeding_tags.push_back("_SP_" + boost::str(boost::format(format)%i));
            }
        }

        if ( seeding_positions->size() >= 4294967296 ) {
            utility_exit_with_message("Too many seeding positions!!!!");
        }

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

		// longxing: as derrick wants to use the IDENTITY director for his C2 ligand, so check this.
		// this is a temporary solution, so the function should not be called when using the IDENTITY director.
		if ( fname == "IDENTITY" ) {
				return true;
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
    bool seeding_by_patchdock,
		float patchdock_min_sasa,
		int patchdock_top_ranks
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

												// remove bad patchdock seeding pos based on the sasa
												if ( utility::string2float(splt[7]) < patchdock_min_sasa ) continue;
												if ( utility::string2int  (splt[1]) > patchdock_top_ranks) continue;

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
        

    if (seeding_positions.size() == 0) {
        utility_exit_with_message( "Error!!!!! Unable to parse seeding file: " + fname );
    }

    return true;
}
    

void dump_xform_file(
    DirectorBase const & director,
    ScenePtr const & scene_minimal,
    double cart_radius,
    double cart_resl,
    double angle_radius,
    double angle_resl
) {

    std::string fname = boost::str(boost::format("xform_pos_cart_rad%.1f_by%.1f_angle%.2f_by%.2f.x")%cart_radius%cart_resl%angle_radius%angle_resl);
    utility::io::ozstream xform_pos( fname );

    std::cout << "Dumping xform_pos to: " << fname << std::endl;

    RifDockIndex rdi;
    uint64_t nest_size = director->size(0, RifDockIndex()).nest_index;

    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for (uint64_t i = 0; i < nest_size; ++i) {
        rdi.nest_index = i;
    
        bool director_success = director->set_scene( rdi, 0, *scene_minimal );
    
        if ( director_success ) {
            EigenXform p = scene_minimal->position(1);
            double ang = Eigen::AngleAxisf( p.rotation() ).angle();
            double ang_degrees = 180.0/3.1415 * ang;

            if ( ang_degrees > angle_radius ) continue;
            #pragma omp critical
            {
            xform_pos << "SP" << " " << i << " " << ang << " "
                      << p.linear().row(0) << " " << p.linear().row(1)<< " " << p.linear().row(2) << " "
                      << p.translation().x() << " " << p.translation().y() << " " << p.translation().z() << std::endl;
          }
        }
    }

    xform_pos.close();

}




}}






