// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/RifFactory.hh>

#include <riflib/util.hh>
#include <riflib/rifdock_typedefs.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/format.hh>
#include <boost/format.hpp>

#include <riflib/rif/RifAccumulators.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <riflib/BurialManager.hh>


#include <scheme/objective/hash/XformMap.hh>
#include <scheme/objective/storage/RotamerScores.hh>

#include <scheme/actor/Atom.hh>
#include <scheme/actor/BackboneActor.hh>
#include <scheme/actor/VoxelActor.hh>
#include <scheme/kinematics/Scene.hh>
#include <scheme/nest/pmap/OriTransMap.hh>
#include <scheme/objective/ObjectiveFunction.hh>
#include <scheme/search/HackPack.hh>

#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <complex>

#include <random>
#include<boost/random/uniform_real.hpp>

#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#endif


namespace devel {
namespace scheme {


template<class C> void call_sort_rotamers( C & v ){ v.second.sort_rotamers(); }
template<class C> void assert_is_sorted  ( C & v ){ runtime_assert( v.second.is_sorted() ); }

template<class XmapIter>
struct XmapKeyIterHelper : public KeyIterHelperBase<typename RifBase::Key> {
    typedef typename RifBase::Key Key;
    XmapKeyIterHelper(XmapIter iter) : iter_(iter) {}
    Key get_key() const override { return iter_->first; }
    void next() override { ++iter_; }
    bool equal(KeyIterHelperBase const & that) const override {
        XmapKeyIterHelper const & concrete = static_cast<XmapKeyIterHelper const &>(that);
        return iter_ == concrete.iter_;
    }
    XmapIter iter_;
};

template< class XMap >
class RifWrapper : public RifBase {

	shared_ptr<XMap> xmap_ptr_;

public:

private:
	RifWrapper() : RifBase() {}
public:
	RifWrapper( shared_ptr<XMap> xmap_ptr, std::string type ) : RifBase(type), xmap_ptr_(xmap_ptr) {}

	virtual bool load( std::istream & in , std::string & description )
	{
		size_t s;
		in.read((char*)&s,sizeof(size_t));
		char buf[9999];
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read(buf,s);
		std::string type_in(buf);
		runtime_assert_msg( type_in == type_, "mismatched rif_types, expected: '" + type_ + "' , got: '" + type_in + "'" );
		return xmap_ptr_->load( in , description );
	}
	virtual bool save( std::ostream & out, std::string & description ) {
		size_t s = type_.size();
		out.write((char*)&s,sizeof(size_t));
		out.write(type_.c_str(),s);
		return xmap_ptr_->save( out, description );
	}

	virtual bool get_xmap_ptr( boost::any * any_p )	{
		bool is_compatible_type =     boost::any_cast< shared_ptr<XMap> const>( any_p );
		if( is_compatible_type ) *any_p = static_cast< shared_ptr<XMap> const>( xmap_ptr_ );
		return is_compatible_type;
	}
	virtual bool get_xmap_const_ptr( boost::any * any_p ) const	{
		bool is_compatible_type =     boost::any_cast< shared_ptr<XMap const> const>( any_p );
		if( is_compatible_type ) *any_p = static_cast< shared_ptr<XMap const> const>( xmap_ptr_ );
		return is_compatible_type;
	}
	virtual bool set_xmap_ptr( boost::any * any_p )	{
		shared_ptr<XMap> const * tmp = boost::any_cast< shared_ptr<XMap> >( any_p );
		if( !tmp ) return false;
		xmap_ptr_ = *tmp;
		return true;
	}


    //Brian
    virtual std::pair< int, int > get_sat1_sat2( EigenXform const & x, int roti ) const override {
    	std::pair< int, int > sat1_sat2( -1, -1);

    	Key k = get_bin_key(x);


        auto const & rs = (*xmap_ptr_)[k];
        int const Nrots = XMap::Value::N;
        for( int i = 0; i < Nrots; ++i ){
            if( rs.empty(i) ) break;
            if (rs.rotamer(i) != roti) continue; 
            std::vector<int> sats;
        	rs.rotamer_sat_groups( i, sats );
        	if ( sats.size() == 0) {
        		return sat1_sat2;
        	}
        	sat1_sat2.first = sats[0];
        	if ( sats.size() == 1 ) {
        		return sat1_sat2;
        	}
        	sat1_sat2.second = sats[1];
        	return sat1_sat2;
        }

        return sat1_sat2;
    }

	size_t size() const override { return xmap_ptr_->size(); }
	float load_factor() const override { return xmap_ptr_->map_.size()*1.f/xmap_ptr_->map_.bucket_count(); }
	size_t mem_use()    const override { return xmap_ptr_->mem_use(); }
	float cart_resl()   const override { return xmap_ptr_->cart_resl_; }
	float ang_resl()    const override { return xmap_ptr_->ang_resl_; }
	size_t sizeof_value_type() const override { return sizeof(typename XMap::Map::value_type); }
	bool  has_sat_data_slots() const override { return XMap::Value::RotScore::UseSat; }

	// will resize to accomodate highest number rotamer
	void get_rotamer_ids_in_use( std::vector<bool> & using_rot ) const override
	{
		typedef typename XMap::Map::value_type MapPair;
		typedef typename MapPair::second_type RotScores;
		for( auto const & v : xmap_ptr_->map_ ){
			RotScores const & xmrot( v.second );
			for( int i = 0; i < RotScores::N; ++i ){
				if( xmrot.empty(i) ) break;
				if( xmrot.rotamer(i) >= using_rot.size() ) using_rot.resize( xmrot.rotamer(i)+1 , false );

				using_rot[ xmrot.rotamer(i) ] = true;
			}
		}

	}

    Key get_bin_key( EigenXform const & x) const override {
        return xmap_ptr_->get_key(x);
    }

    EigenXform get_bin_center( Key const & k) const override {
        return xmap_ptr_->get_center(k);
    }

    void get_rotamers_for_key( Key const & k, std::vector< std::pair< float, int > > & rotscores ) const override{
        typename XMap::Value const & rs = (*xmap_ptr_)[k];
        int const Nrots = XMap::Value::N;
        for( int i = 0; i < Nrots; ++i ){
            if( rs.empty(i) ) break;
            rotscores.push_back( std::make_pair<float,int>( rs.score(i), rs.rotamer(i) ) );
        }
    }

	void
	get_rotamers_for_xform(
		EigenXform const & x,
		std::vector< std::pair< float, int > > & rotscores
	) const override {
        Key const & k = xmap_ptr_->get_key(x);
        get_rotamers_for_key(k, rotscores);
	}

	void finalize_rif() override {
		// sort the rotamers in each cell so best scoring is first
		__gnu_parallel::for_each( xmap_ptr_->map_.begin(), xmap_ptr_->map_.end(), call_sort_rotamers<typename XMap::Map::value_type> );
	}

	// void super_print( std::ostream & out, shared_ptr< RotamerIndex > rot_index_p ) const override { xmap_ptr_->super_print( out, rot_index_p  ); }
	void print( std::ostream & out ) const override { out << (*xmap_ptr_) << std::endl; }
	std::string value_name() const override { return XMap::Value::name(); }

	void collision_analysis( std::ostream & out ) const override {
		using namespace ObjexxFCL::format;
		typedef typename XMap::Value XMapVal;
		int64_t rif_num_collisions  [ XMapVal::N ];
		double  rif_avg_scores      [ XMapVal::N ];
		int64_t rif_avg_scores_count[ XMapVal::N ];
		for( int i = 0; i < XMapVal::N; ++i ){ rif_num_collisions[i]=0; rif_avg_scores[i]=0; rif_avg_scores_count[i]=0; }
		for( auto const & v : xmap_ptr_->map_ ){
			for( int i = 0; i < XMapVal::N; ++i ){
				bool not_empty = !v.second.rotscores_[i].empty();
				if( not_empty ){
					rif_num_collisions[i] += 1;
					rif_avg_scores[i] += v.second.rotscores_[i].score();
					// std::out << v.second.rotscores_[i].score() << std::endl; // WHY SOME WAY TOO LOW?????? fixed.
					rif_avg_scores_count[i]++;
				}
			}
		}
		for( int i = 0; i < XMapVal::N; ++i ) rif_avg_scores[i] /= rif_avg_scores_count[i];

		// out << "======================================================================" << std::endl;
		out << "======================= frac collisions in map =======================" << std::endl;
		out << "======================= cart: "<< xmap_ptr_->cart_resl_ <<
						", ori: " <<  xmap_ptr_->ang_resl_ <<
						" =======================" << std::endl;
		out << "======================================================================" << std::endl;
		float Ecollision = 0.0;
		for( int i = 0; i < XMapVal::N; ++i ){
			float colfrac = rif_num_collisions[i]*1.0/xmap_ptr_->map_.size();
			out << "   Nrots " << I(3,i+1) << " " << F(7,5,colfrac) << " " << F(7,3,rif_avg_scores[i]) << " " << rif_avg_scores_count[i] << std::endl;
			if( i > 0 ){
				float pcolfrac = rif_num_collisions[i-1]*1.0/xmap_ptr_->map_.size();
				Ecollision += i * (pcolfrac-colfrac);
			}
		}
		Ecollision += rif_num_collisions[XMapVal::N-1]*1.0/xmap_ptr_->map_.size() * XMapVal::N;
		out << "E(collisions) = " << Ecollision << std::endl;
		out << "======================================================================" << std::endl;

	}

    RifBaseKeyRange key_range() const override {
        auto b = std::make_shared<XmapKeyIterHelper<typename XMap::Map::const_iterator>>(
            ((typename XMap::Map const &)xmap_ptr_->map_).begin()  );
        auto e = std::make_shared<XmapKeyIterHelper<typename XMap::Map::const_iterator>>(
            ((typename XMap::Map const &)xmap_ptr_->map_).end()  );
        return RifBaseKeyRange(RifBaseKeyIter(b), RifBaseKeyIter(e));
    }

        
        // randomly dump rif residues defined by res_names, and "*" means all 20 amino acids.
        bool random_dump_rotamers( std::vector< std::string > res_names, std::string const file_name, float dump_fraction, shared_ptr<RotamerIndex> rot_index_p ) const override
        {
            std::mt19937 rng(time(0));
            boost::uniform_real<> uniform;
            bool dump_all = false;
            if ( std::find( res_names.begin(), res_names.end(), "*" ) != res_names.end() ) dump_all = true;
            
            const RifBase * base = this;
            shared_ptr<XMap const> from;
            base->get_xmap_const_ptr( from );
            
            
            utility::io::ozstream fout( file_name );
            int64_t count = 1;
            for (auto const & v : from->map_ )
            {
                // this is the position of this grid, will not be used here.
                EigenXform x = from->hasher_.get_center( v.first );
                typename XMap::Value const & rotscores = v.second;
                static int const Nrots = XMap::Value::N;
                for( int i_rs = 0; i_rs < Nrots; ++i_rs )
                {
                    if( rotscores.empty(i_rs) ) {
                        break;
                    }
                    int irot = rotscores.rotamer(i_rs);
                    float score = rotscores.score(i_rs);
                    if ( dump_all || std::find(res_names.begin(), res_names.end(), rot_index_p->rotamers_[irot].resname_) != res_names.end() )
                    {
                        if ( uniform(rng) <= dump_fraction )
                        {
                            fout << std::string("MODEL") << " " << boost::str(boost::format("%.3f")%score) << std::endl;
                            BOOST_FOREACH( SchemeAtom a, rot_index_p->rotamers_.at( irot ).atoms_ ){
                                a.set_position( x * a.position() );
                                a.nonconst_data().resnum = count;
                                a.nonconst_data().chain = 'A';
                                ::scheme::actor::write_pdb( fout, a, nullptr );
                            }
                            fout << std::string("ENDMDL") << std::endl;
                            ++count;
                        }
                    }
                }
            }
            fout.close();
            
        }
        



    // This looks for rifgen rotamers that have their N, CA, CB, and last atom within dump_dist of the residue given
    bool dump_rotamers_near_res( core::conformation::Residue const & res, std::string const & file_name, 
                                        float dump_dist, float dump_frac, shared_ptr<RotamerIndex> rot_index_p ) const override {

        std::string name3 = res.name3();


        std::pair<int, std::string> cb = ( name3 == "GLY" ? std::pair<int, std::string> {1, "CA"} : std::pair<int, std::string>{3, "CB"});
        // rif-number, rosetta number
        std::vector<std::pair<int, std::string>> align_pairs { 
            {0, "N"},
            {1, "CA"},
            cb
        };

        // prepare listed atoms
        std::vector<Eigen::Vector3f> scaff_atoms;
        for ( int i = 0; i < align_pairs.size(); i++ ) {
            std::pair<int, std::string> pair = align_pairs[i];
            numeric::xyzVector<core::Real> xyz = res.xyz(pair.second);
            Eigen::Vector3f vec;
            vec[0] = xyz.x();
            vec[1] = xyz.y();
            vec[2] = xyz.z();
            scaff_atoms.push_back(vec);
        }
        // prepare last atom
        numeric::xyzVector<core::Real> xyz = res.xyz(res.natoms());
        Eigen::Vector3f vec;
        vec[0] = xyz.x();
        vec[1] = xyz.y();
        vec[2] = xyz.z();
        scaff_atoms.push_back(vec);

    	std::cout << "Looking for rotamers within " << dump_dist << "A of in input pdb and of aa " << name3 << std::endl;
    	const RifBase * base = this;
		shared_ptr<XMap const> from;
		base->get_xmap_const_ptr( from );


		float coarse_dist_sq = (dump_dist + 5) * (dump_dist + 5);
		float dump_dist_sq = dump_dist * dump_dist;

		// transform and irot
		std::vector<std::pair<EigenXform, std::pair<int, float>>> to_dump;
		to_dump.reserve( from->map_.size() );


		std::pair<int,int> index_bounds = rot_index_p->index_bounds( name3 );


		for( auto const & v : from->map_ ){
			EigenXform x = from->hasher_.get_center( v.first );

			float dist_sq = (x.translation() - scaff_atoms[2]).squaredNorm();

            if (dist_sq > coarse_dist_sq) continue;

            typename XMap::Value const & rotscores = from->operator[]( x );
            static int const Nrots = XMap::Value::N;
            for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
                if( rotscores.empty(i_rs) ) {
                    break;
                }

                int irot = rotscores.rotamer(i_rs);
                if (irot < index_bounds.first || irot >= index_bounds.second ) continue;

                // check listed atoms
                bool all_good = true;
                for ( int i = 0; i < align_pairs.size(); i++ ) {
                    int rif_atom_no = align_pairs[i].first;
                    Eigen::Vector3f scaff_vec = scaff_atoms[i];


                    SchemeAtom atom = rot_index_p->rotamers_.at( irot ).atoms_[rif_atom_no];
                    Eigen::Vector3f atom_vec = x * atom.position();
                    dist_sq = (atom_vec - scaff_vec).squaredNorm();

                    if (dist_sq > dump_dist_sq) {
                        all_good = false;
                        break;
                    }
                }

                if ( ! all_good ) continue;

                // check last atom
                Eigen::Vector3f scaff_vec = scaff_atoms.back();
                SchemeAtom atom = rot_index_p->rotamers_.at( irot ).atoms_.back();
                Eigen::Vector3f atom_vec = x * atom.position();
                dist_sq = (atom_vec - scaff_vec).squaredNorm();

                if (dist_sq > dump_dist_sq) continue;
                

				float score = rotscores.score(i_rs);

				to_dump.push_back(std::pair<EigenXform, std::pair<int, float>>(x, std::pair<int, float>(irot, score)));
			
			
			}
		}


        if (to_dump.size() == 0) {
            std::cout << "No rotamers found!!!!" << std::endl;
            return false;
        }

		uint64_t num_dump = to_dump.size() * dump_frac;
		uint64_t dump_every = to_dump.size() / num_dump;


		std::cout << "Found " << to_dump.size() << " rotamers. Dumping " << num_dump << " to " << file_name << " ..." << std::endl;


		utility::io::ozstream out( file_name );
		uint64_t dumped = 0;
		for ( uint64_t i = 0; i < to_dump.size(); i ++ ) {
			if ( i % dump_every != 0 ) {
				continue;
			}
			EigenXform x = to_dump[i].first;
			auto inner_pair = to_dump[i].second;
			int irot = inner_pair.first;
			float score = inner_pair.second;


			out << std::string("MODEL") << " " << boost::str(boost::format("%.3f")%score) << std::endl;

            BOOST_FOREACH( SchemeAtom a, rot_index_p->rotamers_.at( irot ).atoms_ ){
                a.set_position( x * a.position() ); 
                a.nonconst_data().resnum = dumped;
                a.nonconst_data().chain = 'A';
                ::scheme::actor::write_pdb( out, a, nullptr );
            }

            dumped ++;

			out << std::string("ENDMDL") << std::endl;

		}

		out.close();


        return true;

    }

    // dumps everything at this bin center that's the same residue type
    void
    dump_rotamers_at_bin_center( 
        core::conformation::Residue const & res,
        std::string const & file_name,
        shared_ptr<RotamerIndex> rot_index_p
    ) const override {

        std::string name3 = res.name3();
        std::pair<int,int> index_bounds = rot_index_p->index_bounds( name3 );

        std::cout << "Dumping residues of type " << name3 << " at bin center to " << file_name << std::endl;

        BBActor bb( res );

        const RifBase * base = this;
        shared_ptr<XMap const> xmap;
        base->get_xmap_const_ptr( xmap );

        EigenXform center = xmap->get_center(xmap->get_key(bb.position()));

        std::vector<std::pair<EigenXform, std::pair<int, float>>> to_dump;

        typename XMap::Value const & rotscores = xmap->operator[]( center );
        static int const Nrots = XMap::Value::N;
        for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
            if( rotscores.empty(i_rs) ) {
                continue;
            }

            int irot = rotscores.rotamer(i_rs);
            if (irot < index_bounds.first || irot >= index_bounds.second ) continue;

            float score = rotscores.score(i_rs);

            to_dump.push_back(std::pair<EigenXform, std::pair<int, float>>(center, std::pair<int, float>(irot, score)));

        }

        if (to_dump.size() == 0) {
            std::cout << "No rotamers found!!!!" << std::endl;
            return;
        }

        std::cout << "Found " << to_dump.size() << " rotamers. Dumping " << to_dump.size() << " to " << file_name << " ..." << std::endl;

        utility::io::ozstream out( file_name );
        uint64_t dumped = 0;
        for ( uint64_t i = 0; i < to_dump.size(); i ++ ) {

            EigenXform x = to_dump[i].first;
            auto inner_pair = to_dump[i].second;
            int irot = inner_pair.first;
            float score = inner_pair.second;


            out << std::string("MODEL") << " " << boost::str(boost::format("%.3f")%score) << std::endl;

            BOOST_FOREACH( SchemeAtom a, rot_index_p->rotamers_.at( irot ).atoms_ ){
                a.set_position( x * a.position() ); 
                a.nonconst_data().resnum = dumped;
                a.nonconst_data().chain = 'A';
                ::scheme::actor::write_pdb( out, a, nullptr );
            }

            dumped ++;

            out << std::string("ENDMDL") << std::endl;

        }

        out.close();


    }

    void
    dump_rifgen_text( EigenXform const & xform, shared_ptr< XMap const > rif, shared_ptr<RotamerIndex> rot_index_p ) const {

        std::cout << xform.translation().transpose() << std::endl;

        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;

        typename XMap::Value const & rotscores = rif->operator[]( xform );
        static int const Nrots = XMap::Value::N;
        for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
            if( rotscores.empty(i_rs) ) {
                std::cout << I(2, i_rs) << std::endl;
                continue;
            }

            int irot = rotscores.rotamer(i_rs);
            std::string oneletter = rot_index_p->oneletter(irot);
            float score = rotscores.score(i_rs);

            std::vector<int> sats;
            rotscores.rotamer_sat_groups(i_rs, sats);
            int sat1 = -1;
            int sat2 = -1;
            if (sats.size() > 0) {
                sat1 = sats[0];
                if (sats.size() > 1) {
                    sat2 = sats[1];
                }
            }


            std::cout << I(2, i_rs) << " " << oneletter << " " << I(3, irot) << " " << F(6,1, score) << " sat1: " << I(4, sat1)
                        << " sat2: " << I(4, sat2) << std::endl;

        }
        std::cout << std::endl;
    }


    // This looks for rifgen bins that are within dump_distance of the res stub
    bool dump_rifgen_text_near_res( core::conformation::Residue const & res, 
                                        float dump_dist, shared_ptr<RotamerIndex> rot_index_p ) const override {

        using ObjexxFCL::format::F;

        // numeric::xyzVector<core::Real> _n  = res.xyz("N");
        // numeric::xyzVector<core::Real> _ca = res.xyz("CA");
        // numeric::xyzVector<core::Real> _c  = res.xyz("C");

        // Eigen::Vector3f n;  n[0]  = _n[0];  n[1]  = _n[1];  n[2] =  _n[2];
        // Eigen::Vector3f ca; ca[0] = _ca[0]; ca[1] = _ca[1]; ca[2] = _ca[2];
        // Eigen::Vector3f c;  c[0]  = _c[0];  c[1]  = _c[1];  c[2] =  _c[2];


        BBActor bb( res );

        std::cout << bb.position().translation().transpose() << std::endl;

        const RifBase * base = this;
        shared_ptr<XMap const> xmap;
        base->get_xmap_const_ptr( xmap );

        std::cout << "Distance 0.00:" << std::endl;

        EigenXform stored_bin = xmap->get_center(xmap->get_key(bb.position()));

        dump_rifgen_text(stored_bin, xmap , rot_index_p );


        EigenXform bbinv = bb.position().inverse();

        int search_points = 0;
        for( auto const & v : xmap->map_ ){
            search_points ++;
            EigenXform x = xmap->hasher_.get_center( v.first );

            float dist = xform_magnitude( bbinv * x, 1 );
            // float dist = (x.translation() - bb.position().translation()).norm();

            if ( dist > dump_dist ) {
                continue;
            }

            std::cout << "Distance: " << F(5, 2, dist) << std::endl;
            dump_rifgen_text( x, xmap , rot_index_p );

        }

        std::cout << "Searched " << search_points << " points" << std::endl;


    }





};



std::string get_rif_type_from_file( std::string fname )
{
	runtime_assert( utility::file::file_exists(fname) );
	utility::io::izstream in( fname );
	runtime_assert( in.good() );
	size_t s;
	in.read((char*)&s,sizeof(size_t));
	char buf[9999];
	for(int i = 0; i < 9999; ++i) buf[i] = 0;
	in.read(buf,s);
	return std::string(buf);
}





//////////////////////////////////// ScoreBBActorVsRIF ////////////////////////////////////////////

	typedef int32_t intRot;



	// should move to libraries somewhere
	// empty class, serves only to position RIF in absolute space
	// struct RIFAnchor {
	//  RIFAnchor() {}
	// };
	std::ostream & operator<<( std::ostream & out, RIFAnchor const& va ){
		return out << "RIFAnchor";
	}
	template< class MetaData >
	void write_pdb( std::ostream & , RIFAnchor const &, MetaData const & ){}

	struct ScoreBBActorvsRIFResult {
		float val_;
		std::vector< std::pair<intRot,intRot> > rotamers_;
		ScoreBBActorvsRIFResult() : val_(0) {}
		ScoreBBActorvsRIFResult( float f ) : val_(f) {}
		operator float() const { return val_; }
		void operator=( float f ) { val_ = f; }
		void operator+=( float f ) { val_ += f; }
		bool operator<( ScoreBBActorvsRIFResult const & other ) const { return val_ < other.val_; }
	};
	struct ScoreBBActorvsRIFScratch {
		shared_ptr< ::scheme::search::HackPack> hackpack_;
		std::vector<bool> is_satisfied_;
		std::vector<bool> has_rifrot_;

        std::vector<bool> requirements_satisfied_;
		std::vector<std::vector<float> > const * rotamer_energies_1b_ = nullptr;
		std::vector< std::pair<int,int> > const * scaffold_rotamers_ = nullptr;
		shared_ptr< BurialManager > burial_manager_;
		shared_ptr< UnsatManager > unsat_manager_;
        shared_ptr< BurialVoxelArray > scaff_burial_grid_;
		shared_ptr<::scheme::objective::storage::TwoBodyTable<float> const> reference_twobody_;
        //std::vector<std::vector<bool>> allowed_irots_;
        shared_ptr<std::vector<std::vector<bool>>> allowed_irots_;
		// sat group vector goes here
		//std::vector<float> is_satisfied_score_;
	};

	template< class BBActor, class RIF, class VoxelArrayPtr >
	struct ScoreBBActorVsRIF
	{
		typedef boost::mpl::true_ HasPre;
		typedef boost::mpl::true_ HasPost;
		typedef ScoreBBActorvsRIFScratch Scratch;
		typedef ScoreBBActorvsRIFResult Result;
		typedef std::pair<RIFAnchor,BBActor> Interaction;
		bool packing_ = false;
		::scheme::search::HackPackOpts packopts_;
		int n_sat_groups_ = 0, require_satisfaction_ = 0, require_n_rifres_ = 0, require_hydrophobic_residue_contacts_ = 0;
        float hydrophobic_ddg_cut_ = 0;
		std::vector< shared_ptr< ::scheme::search::HackPack> > packperthread_;
		std::vector< shared_ptr< BurialManager> > burialperthread_;
		std::vector< shared_ptr< UnsatManager > > unsatperthread_;
        shared_ptr< HydrophobicManager> hydrophobic_manager_;
        
        // the requirements code
        std::vector< int > requirements_;
	private:
		shared_ptr<RIF const> rif_ = nullptr;
	public:
		VoxelArrayPtr target_proximity_test_grid_ = nullptr;
		RifScoreRotamerVsTarget rot_tgt_scorer_;
		std::vector<int> always_available_rotamers_;

		ScoreBBActorVsRIF() {}

		 void clear() {
			packperthread_.clear();
		 }

		static std::string name(){ return "ScoreBBActorVsRIF"; }

		void set_rif( shared_ptr< ::devel::scheme::RifBase const> rif_ptr ){
			rif_ptr->get_xmap_const_ptr( rif_ );
		}

		void init_for_packing(
			// ::scheme::objective::storage::TwoBodyTable<float> const & twob,
			shared_ptr< ::devel::scheme::RotamerIndex > rot_index_p,
            RifScoreRotamerVsTarget const & rot_tgt_scorer,
			::scheme::search::HackPackOpts const & hackpackopts
		){
            rot_tgt_scorer_ = rot_tgt_scorer;
			packing_ = true;
			packperthread_.clear();
			for( int i  = 0; i < ::devel::scheme::omp_max_threads_1(); ++i ){
				shared_ptr< ::scheme::search::HackPack> tmp = make_shared< ::scheme::search::HackPack>(hackpackopts,rot_index_p->ala_rot(),i);
				packperthread_.push_back( tmp );
			}

			packopts_ = hackpackopts;
			always_available_rotamers_.clear();
			for( int irot = 0; irot < rot_index_p->n_primary_rotamers(); ++irot ){
				switch( packopts_.always_available_rotamers_level ){
					case 2:
						always_available_rotamers_.push_back(irot);
						break;
					case 1:
						if( rot_index_p->resname(irot) == "VAL" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "ILE" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "LEU" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "MET" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "PHE" ) always_available_rotamers_.push_back(irot);
						break;
					case 0:
						break;
					default:
						utility_exit_with_message("unknown always_available_rotamers level");
				}
			}
		}
		void init_for_burial(
			shared_ptr< BurialManager > burial_manager,
			shared_ptr< UnsatManager > unsat_manager
		) {
			for( int i  = 0; i < ::devel::scheme::omp_max_threads_1(); ++i ){
				burialperthread_.push_back( burial_manager->clone() );
				unsatperthread_.push_back( unsat_manager->clone() );
			}
		}

		template<class Scene, class Config>
		void pre( Scene const & scene, Result & result, Scratch & scratch, Config const & config ) const
		{

			// Added by brian ////////////////////////
			ScaffoldDataCacheOP data_cache = scene.conformation_ptr(1)->cache_data_;
            runtime_assert( data_cache );
			scratch.rotamer_energies_1b_ = data_cache->local_onebody_p.get();
            scratch.scaff_burial_grid_ = data_cache->burial_grid;
            scratch.allowed_irots_ = data_cache->allowed_irot_at_ires_p;

			//////////////////////////////////////////

			runtime_assert( rif_ );
			runtime_assert( scratch.rotamer_energies_1b_ );
			if( n_sat_groups_ > 0 && burialperthread_.size() == 0 ){
				scratch.is_satisfied_.resize(n_sat_groups_,false); // = new bool[n_sat_groups_];
				for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_[i] = false;
				//scratch.is_satisfied_score_.resize(n_sat_groups_,0.0);
				//for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_score_[i] = 0;
			}
            if ( requirements_.size() > 0 )
            {
                int64_t max = -9e9;
                for( auto val : requirements_ ){ if (val > max ) max = val; }
                scratch.requirements_satisfied_.resize(max+1);
                for( int i = 0; i <= max; ++i ) scratch.requirements_satisfied_[i] = false;
            }
			scratch.has_rifrot_.resize(scratch.rotamer_energies_1b_->size(), false);
			for ( int i = 0; i < scratch.has_rifrot_.size(); i++ ) scratch.has_rifrot_[i] = false;

			if ( burialperthread_.size() > 0 ) {
				scratch.burial_manager_ = burialperthread_.at( ::devel::scheme::omp_thread_num() );
				scratch.burial_manager_->reset();
				scratch.unsat_manager_ = unsatperthread_.at( ::devel::scheme::omp_thread_num() );
				scratch.unsat_manager_->reset();

				scratch.is_satisfied_ = scratch.unsat_manager_->get_presatisfied();
			}

			if( !packing_ ) return;

			// Added by brian ////////////////////////
			scratch.scaffold_rotamers_ = data_cache->local_rotamers_p.get();
			//////////////////////////////////////////

			runtime_assert( scratch.scaffold_rotamers_ );
			runtime_assert( rot_tgt_scorer_.rot_index_p_ );
			runtime_assert( rot_tgt_scorer_.target_field_by_atype_.size() == 22 );
			scratch.hackpack_ = packperthread_.at( ::devel::scheme::omp_thread_num() );

			if ( ! scratch.burial_manager_ ) {
				scratch.hackpack_->reinitialize( data_cache->local_twobody_p );
			} else {
				scratch.hackpack_->reinitialize( data_cache->local_twobody_per_thread.at(::devel::scheme::omp_thread_num()) );
				scratch.reference_twobody_ = data_cache->local_twobody_p;
			}

		}

		template<class Config>
		Result operator()( RIFAnchor const &, BBActor const & bb, Scratch & scratch, Config const& c ) const
		{

			if( target_proximity_test_grid_ && target_proximity_test_grid_->at( bb.position().translation() ) == 0.0 ){
				return 0.0;
			}

			const bool want_sats = scratch.burial_manager_;

			typename RIF::Value const & rotscores = rif_->operator[]( bb.position() );
			static int const Nrots = RIF::Value::N;
			int const ires = bb.index_;
			float bestsc = 0.0;
			for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
				if( rotscores.empty(i_rs) ) {
					break;
				}
				
                int irot = rotscores.rotamer(i_rs);
                if ( scratch.allowed_irots_ && ! scratch.allowed_irots_->at(ires)[irot] ) continue;

				float const rot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
				float score_rot_v_target = rotscores.score(i_rs);

				bool rotamer_satisfies = rotscores.do_i_satisfy_anything(i_rs);

				if( packing_ && packopts_.packing_use_rif_rotamers ){

					if( rot1be <= packopts_.rotamer_onebody_inclusion_threshold || rotamer_satisfies){
						
                        // Very important!!! Do not fill in scratch.is_satisfied_ here!!!
                        //  It needs to collect the BBHbond sats which are unconditionally satisifed
						int sat1 = -1, sat2 = -1, hbcount = 0;
						float const recalc_rot_v_tgt = packopts_.rescore_rots_before_insertion ? 
														rot_tgt_scorer_.score_rotamer_v_target_sat( 
																irot, bb.position(), sat1, sat2, want_sats, hbcount, 10.0, 4 ) :
														score_rot_v_target;

						score_rot_v_target = recalc_rot_v_tgt;
				
						if (( score_rot_v_target + rot1be < packopts_.rotamer_inclusion_threshold &&
							  score_rot_v_target          < packopts_.rotamer_inclusion_threshold ) || rotamer_satisfies){

							float sat_bonus = 0;
							if (rotamer_satisfies) {
                                
                                if ( requirements_.size() > 0 && std::find( requirements_.begin(), requirements_.end(), rotscores.get_requirement_num(i_rs) ) != requirements_.end() ) {
                                    sat_bonus -= 10.0;
                                }
                                
								sat_bonus = packopts_.user_rotamer_bonus_per_chi * rot_tgt_scorer_.rot_index_p_->nchi(irot) +
											packopts_.user_rotamer_bonus_constant;
							}
							if ( ! scratch.burial_manager_ ) scratch.hackpack_->add_tmp_rot( ires, irot, score_rot_v_target + rot1be + sat_bonus );
							else                    scratch.unsat_manager_->add_to_pack_rot( ires, irot, score_rot_v_target + rot1be, sat1, sat2 );
							
						}
					}
					if( packopts_.use_extra_rotamers ){
						auto child_rots = rot_tgt_scorer_.rot_index_p_->child_map_.at(irot);
						for( int crot = child_rots.first; crot < child_rots.second; ++crot ){
							float const crot1be = (*scratch.rotamer_energies_1b_).at(ires).at(crot);
							if( crot1be > packopts_.rotamer_onebody_inclusion_threshold ) continue;

							int sat1 = -1, sat2 = -1, hbcount = 0;
							float const recalc_crot_v_tgt = packopts_.rescore_rots_before_insertion ? 
																rot_tgt_scorer_.score_rotamer_v_target_sat( 
																	crot, bb.position(), sat1, sat2, want_sats, hbcount, 10.0, 4 ) :
																score_rot_v_target; // this is certainly the wrong score

							if( recalc_crot_v_tgt + rot1be < packopts_.rotamer_inclusion_threshold &&
								recalc_crot_v_tgt          < packopts_.rotamer_inclusion_threshold ){
								if ( ! scratch.burial_manager_ ) scratch.hackpack_->add_tmp_rot( ires, crot, recalc_crot_v_tgt + crot1be );
								else                    scratch.unsat_manager_->add_to_pack_rot( ires, crot, recalc_crot_v_tgt + crot1be, sat1, sat2 );
							}
						}
					}
				}

				float const score_rot_tot = score_rot_v_target + rot1be;
				//TaYi change the score_rot_tot cutoff for mark_sat_groups to include not so good rotamers 
				if( n_sat_groups_ > 0 && score_rot_tot < 5.0 ){
					rotscores.mark_sat_groups( i_rs, scratch.is_satisfied_ );
				}
                if ( score_rot_tot < 0.0 ) {
                    // std::cout << "Adding " << ires << std::endl;
                    scratch.has_rifrot_[ires] = true;
                }
								// an arbitrary cutoff value.
								if( requirements_.size() > 0 && score_rot_tot < 2 )
								{
										rotscores.mark_sat_groups( i_rs, scratch.requirements_satisfied_ );
								}

				bestsc = std::min( score_rot_tot , bestsc );
				//}
			}

            // This doesn't respect packopts_.rescore_rots_before_insertion. i.e. this gives garbage at low resolution
            //  Theoretically fixable, but correct rotamers may have been pushed out of the rif
            // // add native scaffold rotamers TODO: this is bugged somehow?
			if( packing_ && packopts_.add_native_scaffold_rots_when_packing ){
				for( int irot = scratch.scaffold_rotamers_->at(ires).first; irot < scratch.scaffold_rotamers_->at(ires).second; ++irot ){
					if( scratch.hackpack_->using_rotamer( ires, irot ) ){
						float const rot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
						int sat1 = -1, sat2 = -1, hbcount = 0;
						float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target_sat( 
															irot, bb.position(), sat1, sat2, want_sats, hbcount, 10.0, 4 );
						float const rot_tot_1b = recalc_rot_v_tgt + rot1be;
						if( rot_tot_1b < -1.0 && recalc_rot_v_tgt < -0.1 ){ // what's this logic??
							if ( ! scratch.burial_manager_ ) scratch.hackpack_->add_tmp_rot( ires, irot, rot_tot_1b );
							else                    scratch.unsat_manager_->add_to_pack_rot( ires, irot, rot_tot_1b, sat1, sat2 );
						}
					}
				}
			}

            // This doesn't respect packopts_.rescore_rots_before_insertion. i.e. this gives garbage at low resolution
            //  Theoretically fixable, but correct rotamers may have been pushed out of the rif
			if( packing_ ){
				for( int irot : always_available_rotamers_ ){
					float const irot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
					if( irot1be > packopts_.rotamer_onebody_inclusion_threshold ) continue;
					const bool want_sats = scratch.burial_manager_;
					int sat1 = -1, sat2 = -1, hbcount = 0;
					float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target_sat( 
																		irot, bb.position(), sat1, sat2, want_sats, hbcount, 10.0, 4 );
					if( recalc_rot_v_tgt + irot1be < packopts_.rotamer_inclusion_threshold &&
						recalc_rot_v_tgt           < packopts_.rotamer_inclusion_threshold ){
						if ( ! scratch.burial_manager_ ) scratch.hackpack_->add_tmp_rot( ires, irot, recalc_rot_v_tgt + irot1be );
						else                    scratch.unsat_manager_->add_to_pack_rot( ires, irot, recalc_rot_v_tgt + irot1be, sat1, sat2 );
					}
				}
			}

			return bestsc;
		}

		template<class Scene, class Config>
		void post( Scene const & scene, Result & result, Scratch & scratch, Config const & config ) const
		{
			if( packing_ ){

				::scheme::search::HackPack & packer( *scratch.hackpack_ );

				float unsat_zerobody = 0;
				if ( scratch.burial_manager_ ) {
                    EigenXform scaffold_xform = scene.position(1);
					unsat_zerobody = scratch.unsat_manager_->prepare_packer( packer, 
                        scratch.burial_manager_->get_burial_weights( scaffold_xform, scratch.scaff_burial_grid_ ),
                        scratch.is_satisfied_ );
				}
				
				result.val_ = packer.pack( result.rotamers_ );
				result.val_ += unsat_zerobody;

				if ( scratch.burial_manager_ ) {scratch.unsat_manager_->fix_packer( packer, scratch.reference_twobody_ );
					// runtime_assert(scratch.reference_twobody_->check_equal(*packer.twob_));
				}
				

                if ( hydrophobic_manager_ ) {
                    std::vector<std::pair<intRot, EigenXform>> irot_and_bbpos;
                    for( int i = 0; i < result.rotamers_.size(); ++i ){
                        BBActor const & bb = scene.template get_actor<BBActor>( 1, result.rotamers_[i].first );
                        int irot = result.rotamers_[i].second;

                        irot_and_bbpos.emplace_back( irot, bb.position() );
                    }
                    std::vector<int> hyd_counts, per_irot_counts;
                    float hydrophobic_ddg = 0;
                    bool pass_better_than = true;
                    int hydrophobic_residue_contacts = hydrophobic_manager_->find_hydrophobic_residue_contacts( irot_and_bbpos, hyd_counts, hydrophobic_ddg,
                                                                                    per_irot_counts, pass_better_than );
                    if ( hydrophobic_residue_contacts < require_hydrophobic_residue_contacts_) {
                        result.val_ = 9e9;
                    }
                    if ( hydrophobic_ddg > hydrophobic_ddg_cut_ ) {
                        result.val_ = 9e9;
                    }
                    if ( ! pass_better_than ) {
                        result.val_ = 9e9;
                    }

                }



                if ( requirements_.size() > 0 )
                {
                    for( int ii = 0; ii < scratch.requirements_satisfied_.size(); ++ii ) scratch.requirements_satisfied_[ii] = false;
                    
                    for( int ii = 0; ii < result.rotamers_.size(); ++ii ){
                        BBActor const & bb = scene.template get_actor<BBActor>( 1, result.rotamers_[ii].first );
                        typename RIF::Value const & rotscores = rif_->operator[]( bb.position() );
                        static int const Nrots = RIF::Value::N;
                        for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
                            if( rotscores.rotamer(i_rs) == result.rotamers_[ii].second ){
                                rotscores.mark_sat_groups( i_rs, scratch.requirements_satisfied_ );
                                break;
                            }
                        }
                    }
                }


				if( n_sat_groups_ > 0 ) for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_[i] = false;
				for( int i = 0; i < result.rotamers_.size(); ++i ){
					BBActor const & bb = scene.template get_actor<BBActor>( 1, result.rotamers_[i].first );
					int sat1=-1, sat2=-1, hbcount=0;
					float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target_sat(
									result.rotamers_[i].second, bb.position(), sat1, sat2, n_sat_groups_ > 0, hbcount, 10.0, 4 );
					// todo: should do extra selection here?
					// if( recalc_rot_v_tgt < -1.0 ){
					//  selected_rotamers.push_back( result.rotamers_[i] );
					// }
					if( n_sat_groups_ > 0 ){
						if( sat1 >= 0 ) scratch.is_satisfied_[ sat1 ] = true;
						if( sat2 >= 0 ) scratch.is_satisfied_[ sat2 ] = true;
					}
				}
				// result.rotamers_ = selected_rotamers;

			} else {

				if ( scratch.burial_manager_ ) {
                    EigenXform scaffold_xform = scene.position(1);
					std::vector<float> burial_weights = scratch.burial_manager_->get_burial_weights( scaffold_xform, scratch.scaff_burial_grid_ );
					result.val_ += scratch.unsat_manager_->calculate_nonpack_score( burial_weights, scratch.is_satisfied_ );
				}


			}


			if( n_sat_groups_ > 0 ){

				int nsat = 0;
				
				for( int i = 0; i < n_sat_groups_; ++i ){
					
					nsat += scratch.is_satisfied_[i];
					//result.val_ += scratch.is_satisfied_score_[i];
				}
				// if (nsat >= 4 ){
				//  #pragma omp critical
				//  {
				//  std::cout << config << "     ";
					
				//  for (int i = 0; i < 10; ++i){
				//      std::cout << " "<< scratch.is_satisfied_[i];

				//  }
				//  std::cout << " " << std::endl;
				//  }
				// }
				// std::cout << "here: " << nsat << std::endl;
				// runtime_assert( 0 <= nsat && nsat <= n_sat_groups_ );

				if( nsat - require_satisfaction_ < 0 ){
					result.val_ = 9e9;
				} //else {
					//result.val_ += -4.0f * nsat;
				//}

				// delete scratch.is_satisfied_;
				// scratch.is_satisfied_.clear();
			}


			// Brian - this is closer to working during packing but still doesn't work
			if ( require_n_rifres_ > 0 ) {
				if ( packing_ ) {
                    int num_rots = 0;
					for( int i = 0; i < result.rotamers_.size(); ++i ){
                        int irot = result.rotamers_[i].second;
                        if ( irot > 0 ) {
                            num_rots++;
                        }
					}
					// std::cout << result.rotamers_.size() << " ";
					if ( num_rots < require_n_rifres_ ) {
						result.val_ = 9e9;
					}
				} else {
					int count = 0;
					for ( int i = 0; i < scratch.has_rifrot_.size(); i++ ) {
						if ( scratch.has_rifrot_[i] ) {
							count ++;
						}
					}
					// std::cout << "Found " << count << std::endl;
					if (count < require_n_rifres_ ) {
						result.val_ = 9e9;
					}
				}
			}
            
            if ( requirements_.size() > 0 )
            {
                bool pass = true;
                for ( auto const & x : requirements_ ) {
										pass &= scratch.requirements_satisfied_[x];
								}
                if ( !pass ) result.val_ = 9e9;
            }

				// #ifdef USE_OPENMP
				// #pragma omp critical
				// #endif
				// std::cout << result.rotamers_.size() << " " << result.val_ << std::endl;
			// if( result.rotamers_.size() == 0 ){
			//  #ifdef USE_OPENMP
			//  #pragma omp critical
			//  #endif
			//  std::cout << "no rotamers!" << std::endl;
			// }

			// std::cout << "bound score: " << result.val_ << ", packscore: " << packscore << std::endl;
			// for( int i = 0; i < result.rotamers_.size(); ++i ){
			//  std::cout << "res: " << result.rotamers_[i].first << ", rotamer: " << result.rotamers_[i].second << std::endl;
			// }

		}
	};
	template< class B, class X, class V >
	std::ostream & operator<<( std::ostream & out, ScoreBBActorVsRIF<B,X,V> const& si ){ return out << si.name(); }




//////////////////////////////////// ScoreBBHBondActorVsRIF ////////////////////////////////////////////


    struct ScoreBBHBondActorvsRIFResult {
        float val_;
        ScoreBBHBondActorvsRIFResult() : val_(0) {}
        ScoreBBHBondActorvsRIFResult( float f ) : val_(f) {}
        operator float() const { return val_; }
        void operator=( float f ) { val_ = f; }
        void operator+=( float f ) { val_ += f; }
        bool operator<( ScoreBBHBondActorvsRIFResult const & other ) const { return val_ < other.val_; }
    };

    template< class BBHBondActor, class VoxelArrayPtr >
    struct ScoreBBHBondActorVsRIF
    {

        typedef ScoreBBActorvsRIFScratch Scratch;       // IMPORTANT!! This is the same as ScoreBBActorVsRIF
        typedef ScoreBBHBondActorvsRIFResult Result;
        typedef std::pair<RIFAnchor,BBHBondActor> Interaction;
        
        VoxelArrayPtr target_proximity_test_grid_ = nullptr;
        RifScoreRotamerVsTarget rot_tgt_scorer_;
        bool initialized_ = false;
        bool packing_ = false;

        ScoreBBHBondActorVsRIF() {}

         void clear() {
         }

        static std::string name(){ return "ScoreBBHBondActorVsRIF"; }

        void init(
            RifScoreRotamerVsTarget const & rot_tgt_scorer,
            float scaff_bb_hbond_weight,
            bool packing
        ){
            rot_tgt_scorer_ = rot_tgt_scorer;   // This must make a copy!!!!!!
            rot_tgt_scorer_.hbond_weight_ = scaff_bb_hbond_weight;
            initialized_ = true;
            packing_ = packing;
        }

        template<class Config>
        Result operator()( RIFAnchor const &, BBHBondActor const & bbh, Scratch & scratch, Config const& c ) const
        {

            if ( ! initialized_ ) return 0.0; // this is to block lower resolutions

            // if( target_proximity_test_grid_ && target_proximity_test_grid_->at( bbh.hbond_rays().front().horb_cen ) == 0.0 ){
            //     return 0.0;
            // }

            float score = 0;

            if ( true ) {
                int sat1 = -1, sat2 = -1, hbcount = 0;

                if ( bbh.is_donor() ) {
                    score += rot_tgt_scorer_.score_donor_rays_v_target( bbh.hbond_rays(), sat1, sat2, hbcount );
                } else {
                    score += rot_tgt_scorer_.score_acceptor_rays_v_target( bbh.hbond_rays(), sat1, sat2, hbcount );
                }

                if ( scratch.is_satisfied_.size() > 0 ) {
                    if ( sat1 > -1 ) scratch.is_satisfied_.at(sat1) = true;
                    if ( sat2 > -1 ) scratch.is_satisfied_.at(sat2) = true;
                }
            } else {

                typedef typename DonorAcceptorCache::Sat Sat;
                bool any = false;

                // When in doubt on which is which here, just check the callpath in rot_tgt_scorer_
                shared_ptr<DonorAcceptorCache> const & cache = bbh.is_donor() ? rot_tgt_scorer_.target_acceptor_cache_ : rot_tgt_scorer_.target_donor_cache_;
                Sat adder = bbh.is_donor() ? rot_tgt_scorer_.target_donors_.size() : 0  ;

                runtime_assert(cache);
                Eigen::Vector3f super_far_away(1e5, 1e5, 1e5);

                for ( auto const & ray : bbh.hbond_rays() ) {
                    std::vector<Sat>::const_iterator sats_iter =  rot_tgt_scorer_.target_acceptor_cache_->at( ray.horb_cen );

                    Sat don_or_acc = 0;
                    while ( (don_or_acc = *(sats_iter++)) != DonorAcceptorCache::CACHE_MAX_SAT ) {
  
                        don_or_acc += adder;
                        scratch.is_satisfied_.at(don_or_acc ) = true;
                        any = true;
                    }

                } 
                score = any ? -rot_tgt_scorer_.hbond_weight_ : 0;

            }

            return score;
        }

    };
    template< class B, class V >
    std::ostream & operator<<( std::ostream & out, ScoreBBHBondActorVsRIF<B,V> const& si ){ return out << si.name(); }




//////////////////////////////////// ScoreBBSasaActorVsRIF ////////////////////////////////////////////


    struct ScoreBBSasaActorvsRIFResult {
        float val_;
        ScoreBBSasaActorvsRIFResult() : val_(0) {}
        ScoreBBSasaActorvsRIFResult( float f ) : val_(f) {}
        operator float() const { return val_; }
        void operator=( float f ) { val_ = f; }
        void operator+=( float f ) { val_ += f; }
        bool operator<( ScoreBBSasaActorvsRIFResult const & other ) const { return val_ < other.val_; }
    };

    template< class BBSasaActor, class VoxelArrayPtr >
    struct ScoreBBSasaActorVsRIF
    {

        // typedef ScoreBBSasavsRIFScratch Scratch;      
        typedef ScoreBBSasaActorvsRIFResult Result;
        typedef std::pair<RIFAnchor,BBSasaActor> Interaction;

        shared_ptr<BurialVoxelArray> sasa_grid_;
        float threshold_;
        float multiplier_;
        bool initialized_ = false;

        ScoreBBSasaActorVsRIF() {}

         void clear() {
         }

        static std::string name(){ return "ScoreBBSasaActorVsRIF"; }

        void init(
            shared_ptr<BurialVoxelArray> const & sasa_grid,
            float threshold,
            float multiplier
        ){
            sasa_grid_ = sasa_grid;
            threshold_ = threshold;
            multiplier_ = multiplier;
            initialized_ = true;
        }

        template<class Config>
        Result operator()( RIFAnchor const &, BBSasaActor const & bbs, Config const& c ) const
        {
            if ( ! initialized_ ) return 0;     // this is to block lower resolutions

            float score = 0;

            for ( Eigen::Vector3f const & pt : bbs.sasa_points() ) {

                float grid_score = sasa_grid_->at( pt );

                if ( grid_score > threshold_ ) {
                    score += 1;
                }

            }

            score *= multiplier_;

            return score;
        }

    };
    template< class B, class V >
    std::ostream & operator<<( std::ostream & out, ScoreBBSasaActorVsRIF<B,V> const& si ){ return out << si.name(); }




template< class XMap >
struct RifFactoryImpl :
	public RifFactory,
	enable_shared_from_this< RifFactoryImpl<XMap> >
 {


	typedef ScoreBBActorVsRIF<
			BBActor,
			// debugXMap,
			XMap, // no worky worky
			VoxelArrayPtr
		> MyScoreBBActorRIF;

    typedef ScoreBBHBondActorVsRIF<
            BBHBondActor,
            VoxelArrayPtr
        > MyScoreBBHBondActorRIF;

    typedef ScoreBBSasaActorVsRIF<
            BBSasaActor,
            VoxelArrayPtr
        > MyScoreBBSasaActorRIF;


	typedef ::scheme::objective::ObjectiveFunction<
			boost::mpl::vector<          // Do not change this order, only append
				MyScoreBBActorRIF,
				MyClashScore,
                MyScoreBBHBondActorRIF,
                MyScoreBBSasaActorRIF   
			>,
			int // Config type, just resl
		> MyRIFObjective;

	typedef ::scheme::objective::integration::SceneObjectiveParametric<
			ParametricScene,
			MyRIFObjective,
			MyScoreBBActorRIF
		> MySceneObjectiveRIF;

	RifFactoryImpl( RifFactoryConfig const & config ) : RifFactory(config) {}

	virtual RifPtr
	create_rif( float cart_resl=0, float ang_resl=0, float cart_bound=0 ) const
	{
	   if( cart_resl != 0 && ang_resl != 0 && cart_bound != 0 ){
	       return make_shared<RifWrapper<XMap> >( make_shared<XMap>( cart_resl, ang_resl, cart_bound ), this->config().rif_type );
	   } else if( cart_resl == 0 && ang_resl == 0 && cart_bound == 0 ){
	       return make_shared<RifWrapper<XMap> >( make_shared<XMap>(), this->config().rif_type );
	   } else {
	   		utility_exit_with_message("some XformMap constructor values specified, others not!");
	   }
	}

	virtual RifPtr
	create_rif_from_rif( RifConstPtr refrif, float cart_resl, float ang_resl, float cart_bound ) const {
		runtime_assert( this->config().rif_type == refrif->type() );
		RifPtr rif = create_rif( cart_resl, ang_resl, cart_bound );

		shared_ptr<XMap> to;
		rif->get_xmap_ptr( to );
		shared_ptr<XMap const> from;
		refrif->get_xmap_const_ptr( from );

		// std::cout << "create rif progress "; std::cout.flush();


		// old
		int progress0 = 0;
		for( auto const & v : from->map_ ){
			// if( ++progress0 % std::max((size_t)1,(from->size()/100)) == 0 ){
				// std::cout << '*'; std::cout.flush();
			// }
			EigenXform x = from->hasher_.get_center( v.first );

			uint64_t k = to->hasher_.get_key(x);
			typename XMap::Map::iterator iter = to->map_.find(k);
			if( iter == to->map_.end() ){
				to->map_.insert( std::make_pair(k,v.second) );
			} else {
				iter->second.merge( v.second );

			}
		}
		// // std::cout << std::endl;

		// new
		// int progress0 = 0;
		// for( auto const & v : from->map_ ){
		// 	// if( ++progress0 % std::max((size_t)1,(from->size()/100)) == 0 ){
		// 		// std::cout << '*'; std::cout.flush();
		// 	// }
		// 	EigenXform x = from->hasher_.get_center( v.first );
		// 	std::vector<uint64_t> keys = to->hasher_.get_key_and_nbrs(x);
		// 	for ( uint64_t const & k : keys ) {
		// 		typename XMap::Map::iterator iter = to->map_.find(k);
		// 		if( iter == to->map_.end() ){
		// 			to->map_.insert( std::make_pair(k,v.second) );
		// 		} else {
		// 			iter->second.merge( v.second );
		// 		}
		// 	}
		// }
		// std::cout << std::endl;


		return rif;
	}

	virtual	shared_ptr<rif::RifAccumulator>
	create_rif_accumulator( float cart_resl, float ang_resl, float cart_bound, size_t scratchM ) const {
		return make_shared< rif::RIFAccumulatorMapThreaded<XMap> >(
			this->shared_from_this(),
			cart_resl, ang_resl, cart_bound,
			scratchM
		);
	}

	virtual RifPtr
	create_rif_from_file( std::string const & fname, std::string & description ) const
	{
		RifPtr rif = this->create_rif();
		runtime_assert( rif );
		if( ! utility::file::file_exists(fname) ){
			utility_exit_with_message("create_rif_from_file missing file: " + fname );
		}
		utility::io::izstream in( fname );
		if( !in.good() ) return nullptr;
		bool success = rif->load( in, description );
		in.close();
		if( success ) return rif;
		else return nullptr;
	}

	virtual ScenePtr
	create_scene() const {
		ScenePtr s = make_shared<ParametricScene>(2);
		s->add_actor( 0, RIFAnchor() );
		return s;
	}

	virtual bool
	create_objectives(
		RifSceneObjectiveConfig const & config,
		std::vector<ObjectivePtr> & objectives,
		std::vector<ObjectivePtr> & packing_objectives
	) const {

		for( int i_so = 0; i_so < config.rif_ptrs.size(); ++i_so ){
			if( i_so <= config.rif_ptrs.size()-1 ){
				if ( config.rif_ptrs[i_so] == nullptr ) continue;
				shared_ptr< MySceneObjectiveRIF> objective = make_shared<MySceneObjectiveRIF>();
				objective->objective.template get_objective< MyScoreBBActorRIF >().set_rif( config.rif_ptrs[i_so] );
				if( i_so > 0 ){
					int i_tptg = std::min( 2, i_so-1 );
					// std::cout << "resl " << config.resolutions[i_so] << " using target_prox_grid " << config.resolutions[i_tptg] << std::endl;
					// use vdw grids as tgt proximity measure, use CH3 atom
					objective->objective.template get_objective<MyScoreBBActorRIF>().target_proximity_test_grid_ =
						config.target_bounding_by_atype->at(i_tptg).at(5);
				}
				objective->config = i_so;
				objectives.push_back( objective );
			}
		}
		// shared_ptr< MySceneObjectiveRIF> objective = make_shared<MySceneObjectiveRIF>();
		// objective->objective.template get_objective<MyScoreBBActorRIF>().set_rif( config.rif_ptrs.back() );
		// objective->config = config.rif_ptrs.size()-1;
		// objectives.push_back( objective );

        for( int i_so = 0; i_so < config.rif_ptrs.size(); ++i_so ){
				if ( config.rif_ptrs[i_so] == nullptr ) continue;
    		shared_ptr< MySceneObjectiveRIF> packing_objective = make_shared<MySceneObjectiveRIF>();
    		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().set_rif( config.rif_ptrs[i_so] );
    		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).config = i_so;
    		if( config.require_satisfaction > 0 ){
    			dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template
    				get_objective<MyScoreBBActorRIF>().n_sat_groups_ = config.n_sat_groups;
    			dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template
    				get_objective<MyScoreBBActorRIF>().require_satisfaction_ = config.require_satisfaction;
    		}
				if( config.requirements.size() > 0 ){
						dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().requirements_ = config.requirements;
				}
        packing_objectives.push_back( packing_objective );
        }

		for( auto op : objectives ){
			// dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().rotamer_energies_1b_ = config.local_onebody;
			// dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().scaffold_rotamers_ = config.local_rotamers;
			if( config.require_satisfaction > 0 ){
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().n_sat_groups_ = config.n_sat_groups;
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().require_satisfaction_ = config.require_satisfaction;
			}
			if (config.require_n_rifres > 0 ) {
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().require_n_rifres_ = config.require_n_rifres;
			}
			if ( config.requirements.size() > 0 ){
			  dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().requirements_ = config.requirements;
			}
		}
		// dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().rotamer_energies_1b_ = config.local_onebody;
		// dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().scaffold_rotamers_ = config.local_rotamers;

		// use 4.0A vdw grid for CH3 atoms as proximity test

        for( int i_so = 0; i_so < packing_objectives.size(); ++i_so ){

            shared_ptr< MySceneObjectiveRIF> packing_objective = std::dynamic_pointer_cast<MySceneObjectiveRIF>(packing_objectives[i_so]);

            ::scheme::search::HackPackOpts local_packopts = *config.packopts;
            // Only rescore on the full resolution
            local_packopts.rescore_rots_before_insertion = i_so == packing_objectives.size() - 1;

    		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>()
    							.target_proximity_test_grid_ = config.target_bounding_by_atype->at(2).at(5);
    		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().init_for_packing(
    			// *config.local_twobody,
    			config.rot_index_p,
                config.rot_tgt_scorer,
    			local_packopts
    		);
        }

        // Only want to do BBHbonds at highest resl otherwise they are meaningless
        if ( objectives.size() ) {
            dynamic_cast<MySceneObjectiveRIF&>(*objectives.back()).objective.template get_objective<MyScoreBBHBondActorRIF>()
                                    .target_proximity_test_grid_ = config.target_bounding_by_atype->at(2).at(5);
            dynamic_cast<MySceneObjectiveRIF&>(*objectives.back()).objective.template get_objective<MyScoreBBHBondActorRIF>()
                                    .init( config.rot_tgt_scorer, config.scaff_bb_hbond_weight, false );
        }
        if ( packing_objectives.back() ) {
            dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBHBondActorRIF>()
                                    .target_proximity_test_grid_ = config.target_bounding_by_atype->at(2).at(5);
            dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBHBondActorRIF>()
                                    .init( config.rot_tgt_scorer, config.scaff_bb_hbond_weight, true );
        }

        // Only want to do BBSasa at highest resl otherwise they are meaningless
        if ( objectives.size() ) {
            dynamic_cast<MySceneObjectiveRIF&>(*objectives.back()).objective.template get_objective<MyScoreBBSasaActorRIF>()
                                    .init( config.sasa_grid, config.sasa_threshold, config.sasa_multiplier );
        }
        if ( packing_objectives.back() ) {
            dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBSasaActorRIF>()
                                    .init( config.sasa_grid, config.sasa_threshold, config.sasa_multiplier );
        }



        if ( config.burial_manager ) {
            if ( objectives.size() ) {
                dynamic_cast<MySceneObjectiveRIF&>(*objectives.back()).objective.template get_objective<MyScoreBBActorRIF>()
                                .init_for_burial( config.burial_manager, config.unsat_manager );
            }
            if ( packing_objectives.size() ) {
                dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBActorRIF>()
                                .init_for_burial( config.burial_manager, config.unsat_manager );
            }
        }

        if ( config.hydrophobic_manager ) {
            // Only the packing objectives need this
            if ( packing_objectives.size() ) {
                dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBActorRIF>()
                                .hydrophobic_manager_ = config.hydrophobic_manager;
                dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBActorRIF>()
                                .require_hydrophobic_residue_contacts_ = config.require_hydrophobic_residue_contacts;
                dynamic_cast<MySceneObjectiveRIF&>(*packing_objectives.back()).objective.template get_objective<MyScoreBBActorRIF>()
                                .hydrophobic_ddg_cut_ = config.hydrophobic_ddg_cut;
            }
        }

		return true;

	}

// DELETE THIS, ONLY USED WHEN TRAINING SASA
    void
    set_sasa_params(
        std::vector<ObjectivePtr> & objectives,
        shared_ptr<BurialVoxelArray> sasa_grid,
        float sasa_threshold,
        float sasa_multiplier
    ) const override {
        if ( objectives.size() ) {
            dynamic_cast<MySceneObjectiveRIF&>(*objectives.back()).objective.template get_objective<MyScoreBBSasaActorRIF>()
                                    .init( sasa_grid, sasa_threshold, sasa_multiplier );
        }
    }


    std::vector<shared_ptr<UnsatManager>> &
    get_unsatperthread( ObjectivePtr & objective ) const override {
        return dynamic_cast<MySceneObjectiveRIF&>(*objective).objective.template get_objective<MyScoreBBActorRIF>().unsatperthread_;
    }


};


shared_ptr<RifFactory>
create_rif_factory( RifFactoryConfig const & config )
{
	if( config.rif_type == "RotScore" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 12, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 24 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 32 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScore64" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 28, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScore128" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 60, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 120 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 128 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScoreSat" )
	{
		typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 9, -4> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 14, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
    else if( config.rif_type == "RotScoreSat96" )
    {
        typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 9, -4> crfRotScore;
        typedef ::scheme::objective::storage::RotamerScores< 22, crfRotScore > crfXMapValue;
        BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 88 );
        typedef ::scheme::objective::hash::XformMap<
                EigenXform,
                crfXMapValue,
                ::scheme::objective::hash::XformHash_bt24_BCC6
            > crfXMap;
        BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 96 );

        return make_shared< RifFactoryImpl<crfXMap> >( config );
    }
    else if( config.rif_type == "RotScoreSat128" )
    {
        typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 9, -4> crfRotScore;
        typedef ::scheme::objective::storage::RotamerScores< 30, crfRotScore > crfXMapValue;
        BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 120 );
        typedef ::scheme::objective::hash::XformMap<
                EigenXform,
                crfXMapValue,
                ::scheme::objective::hash::XformHash_bt24_BCC6
            > crfXMap;
        BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 128 );

        return make_shared< RifFactoryImpl<crfXMap> >( config );
    }
    else if ( config.rif_type == "RotScoreReq" )
    {
        using SatDatum = ::scheme::objective::storage::SatisfactionDatum<uint8_t>;
        typedef ::scheme::objective::storage::RotamerScoreSat<
        uint16_t, 9, -13, SatDatum, 1> crfRotScore;
        typedef ::scheme::objective::storage::RotamerScores< 28, crfRotScore > crfXMapValue;
        BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 84 );
        typedef ::scheme::objective::hash::XformMap<
        EigenXform,
        crfXMapValue,
        ::scheme::objective::hash::XformHash_bt24_BCC6
        > crfXMap;
        BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 96 );
        
        return make_shared< RifFactoryImpl<crfXMap> >( config );
    }
    else if( config.rif_type == "RotScoreSat_2x16" )
    {
        using SatDatum = ::scheme::objective::storage::SatisfactionDatum<uint16_t>;
        typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 9, -13, SatDatum> crfRotScore;
        typedef ::scheme::objective::storage::RotamerScores< 19, crfRotScore > crfXMapValue;
        // PRINT_SIZE_AS_ERROR<sizeof(crfXMapValue)>()();
        BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 114 );
        typedef ::scheme::objective::hash::XformMap<
                EigenXform,
                crfXMapValue,
                ::scheme::objective::hash::XformHash_bt24_BCC6
            > crfXMap;

        // PRINT_SIZE_AS_ERROR<sizeof(crfXMap::Map::value_type)>()();
        BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 128 );

        return make_shared< RifFactoryImpl<crfXMap> >( config );
    }
	else if( config.rif_type == "RotScoreSat_1x16" )
	{
		using SatDatum = ::scheme::objective::storage::SatisfactionDatum<uint16_t>;
		typedef ::scheme::objective::storage::RotamerScoreSat<
					uint16_t, 9, -13, SatDatum, 1> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 14, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );

	} 
	else if( config.rif_type == "Rot10Score6Sat16" )
	{
		typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 10, -4> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 14, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	} else
	{
		utility_exit_with_message( "create_rif_factory_inner: unknown rif type "+config.rif_type );
	}

}







}}
