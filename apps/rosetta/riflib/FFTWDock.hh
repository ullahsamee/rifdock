// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#ifndef INCLUDED_pilot_will_FFTWDock_hh
#define INCLUDED_pilot_will_FFTWDock_hh



#include <iostream>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <boost/foreach.hpp>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <complex>
#include <utility/io/ozstream.hh>

#include <fftw3.h>


//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// FFTW Dock //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////// welzl ////////////////////
	class Sphere {
	public:
		numeric::xyzVector<core::Real> center;
		core::Real radius;

		Sphere() : center(0,0,0),radius(1) {}
		Sphere(numeric::xyzVector<core::Real> const & c, core::Real r) : center(c),radius(r) {}

		Sphere(numeric::xyzVector<core::Real> const & O);   // numeric::xyzVector<core::Real>-Sphere
		Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A);   // Sphere through two points
		Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A, numeric::xyzVector<core::Real> const & B);   // Sphere through three points
		Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A,
		       numeric::xyzVector<core::Real> const & B, numeric::xyzVector<core::Real> const & C);   // Sphere through four points

		// Distance from p to boundary of the Sphere
		core::Real signdis (numeric::xyzVector<core::Real> const & P) const { return center.distance(P) - radius; }
		core::Real signdis2(numeric::xyzVector<core::Real> const & P) const { return center.distance_squared(P)-radius*radius; } // NOT square of signdis!

		bool intersect( Sphere const & other ) const {
			return center.distance_squared(other.center) <= (radius+other.radius)*(radius+other.radius);
		}
		bool contact( Sphere const & other, core::Real const & contact_dis ) const {
			return center.distance_squared(other.center) <= (radius+other.radius+contact_dis)*(radius+other.radius+contact_dis);
		}
		bool contains( numeric::xyzVector<core::Real> const & pt ) const { return center.distance_squared(pt) < radius*radius; }

		// void dump_pdb(std::ostream & out, io::DumpPdbOpts const & opts, io::DumpPdbState & state ) const;
	};
	// std::ostream & operator<<(std::ostream & out, Sphere const & s);

	Sphere::Sphere(numeric::xyzVector<core::Real> const & O)	{
		radius = 0 + 1e-8;
		center = O;
	}

	Sphere::Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A){
		numeric::xyzVector<core::Real> const a = A - O;
		numeric::xyzVector<core::Real> const o = 0.5 * a;
		radius = o.length() + 1e-8;
		center = O + o;
	}

	Sphere::Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A, numeric::xyzVector<core::Real> const & B){
		numeric::xyzVector<core::Real>  const a = A-O, b = B-O;
		core::Real const det_2 = 2.0 * ((a .cross(b)).dot(a .cross(b)) );
		numeric::xyzVector<core::Real>  const o = (b.dot(b) * ((a.cross(b)).cross(a)) +
		                a.dot(a) * (b.cross(a.cross(b)))) / det_2;
		radius = o.length() + 1e-8;
		center = O + o;
	}

	Sphere::Sphere(numeric::xyzVector<core::Real> const & O, numeric::xyzVector<core::Real> const & A,
	               numeric::xyzVector<core::Real> const & B, numeric::xyzVector<core::Real> const & C){
		numeric::xyzVector<core::Real>  const a = A-O, b = B-O, c = C-O;
		core::Real const det_2 = 2.0 * numeric::xyzMatrix<core::Real>(numeric::xyzMatrix<core::Real>::cols(a,b,c)).det();
		numeric::xyzVector<core::Real>  const o = (c.dot(c) * a.cross(b) +
		                b.dot(b) * c.cross(a) +
		                a.dot(a) * b.cross(c) ) / det_2;
		radius = o.length() + 1e-8;
		center = O + o;
	}

	Sphere welzl_sphere(
		std::vector<numeric::xyzVector<core::Real> > const & points,
		core::Size index,
		std::vector<numeric::xyzVector<core::Real> > & sos,
		core::Size numsos
	){
		// if no input points, the recursion has bottomed out. Now compute an
		// exact sphere based on points in set of support (zero through four points)
		if (index == 0) {
			switch( numsos ){
				case 0: return Sphere();
				case 1: return Sphere(sos[0]);
				case 2: return Sphere(sos[0], sos[1]);
				case 3: return Sphere(sos[0], sos[1], sos[2]);
				case 4: return Sphere(sos[0], sos[1], sos[2], sos[3]);
			}
		}
		// Pick a point at "random" (here just the last point of the input set)
		--index;
		// Recursively compute the smallest bounding sphere of the remaining points
		Sphere smallestSphere = welzl_sphere(points, index, sos, numsos ); // (*)
		// If the selected point lies inside this sphere, it is indeed the smallest
		if( smallestSphere.contains( points[index] ) ) return smallestSphere;
		// Otherwise, update set of support to additionally contain the new point
		runtime_assert( numsos < 4 );
		sos[numsos] = points[index];
		// Recursively compute the smallest sphere of remaining points with new s.o.s.
		return welzl_sphere( points, index, sos, numsos+1 );
	}

	Sphere welzl_sphere(std::vector<numeric::xyzVector<core::Real> > const & points){
		std::vector<numeric::xyzVector<core::Real> > sos(4);
		return welzl_sphere(points,points.size(),sos,0);
	}

	Sphere approx_bounding_sphere(
		std::vector<numeric::xyzVector<core::Real> > const & points
	){
		numeric::xyzVector<core::Real> center(0,0,0);
		core::Real radius = -1;
		if(points.size() > 0)	{

			center = 0;
			for(core::Size i = 0; i < points.size(); i++) center += points[i];
			center /= (core::Real)points.size();

			for(core::Size i = 0; i < points.size(); i++){
				core::Real d2 = points[i].distance_squared(center);
				if(d2 > radius) radius = d2;
			}
			radius = sqrt(radius) + 1e-8;
		}
		return Sphere(center, radius);
	}

struct FFTHit {
	float score;
	numeric::xyzVector<float> trans;
	numeric::xyzMatrix<float> rot1, rot2;
 };
 typedef utility::vector1<FFTHit> FFTHits;
 bool cmpscore( FFTHit const & i, FFTHit const & j) { return i.score > j.score ; }

struct FFTWDock {
	typedef utility::vector1<core::Size> Sizes;
	typedef utility::vector1<std::string> Strings;
	typedef utility::vector1<core::id::AtomID> AIDs;
	typedef ObjexxFCL::FArray3D< std::complex<double> > Complex3D;

	core::pose::Pose pose1_, pose2_;
	numeric::xyzVector<core::Real> center1_, center2_;
	double clash_range_, inter_range_, grid_spacing_, min_contact_volume_;
	AIDs      aids_clash_1_, aids_clash_2_, aids_inter_1_, aids_inter_2_;
	Complex3D grid_clash_1_, grid_clash_2_, grid_inter_1_, grid_inter_2_;
	Complex3D freq_clash_1_, freq_clash_2_, freq_inter_1_, freq_inter_2_;
	fftw_plan fftw_plan_forward_, fftw_plan_reverse_;
	numeric::xyzVector<core::Real> filter_direction_;
	core::Real filter_distance_;
	bool verbose_, plans_made_;
	numeric::xyzMatrix<core::Real> pre_rotation_;
	core::Size gsize_x_,gsize_y_,gsize_z_;
	double radius1, radius2;

	FFTWDock(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2
	){
		pose1_ = pose1;
		pose2_ = pose2;
		pre_rotation_ = numeric::xyzMatrix<core::Real>::identity();
		verbose_ = false;
		plans_made_ = false;
	}

	virtual ~FFTWDock(){
		if(plans_made_){
			fftw_destroy_plan(fftw_plan_forward_);
			fftw_destroy_plan(fftw_plan_reverse_);
		}
	}

	core::Size grid_size_x() const { return gsize_x_; }
	core::Size grid_size_y() const { return gsize_y_; }
	core::Size grid_size_z() const { return gsize_z_; }

	void init(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		double clash_range,
		double inter_range,
		double grid_spacing,
		double min_contact_volume,
		AIDs const & aids_clash_1,
		AIDs const & aids_inter_1,
		AIDs const & aids_clash_2,
		AIDs const & aids_inter_2,
		numeric::xyzVector<core::Real> filter_direction = numeric::xyzVector<core::Real>(0,0,1),
		core::Real filter_distance = 9e9,
		bool fixed = false
	){
		// copy pose coords
		for(core::Size ir = 1; ir <= pose1_.n_residue(); ++ir)
			for(core::Size ia = 1; ia <= pose1_.residue_type(ir).natoms(); ++ia)
				pose1_.set_xyz( core::id::AtomID(ia,ir), pose1.xyz(core::id::AtomID(ia,ir)) );
		for(core::Size ir = 1; ir <= pose2_.n_residue(); ++ir)
			for(core::Size ia = 1; ia <= pose2_.residue_type(ir).natoms(); ++ia)
				pose2_.set_xyz( core::id::AtomID(ia,ir), pose2.xyz(core::id::AtomID(ia,ir)) );

		runtime_assert_msg( pre_rotation_ == numeric::xyzMatrix<core::Real>::identity(), "pre_rotation not implemented (correctly)" );
		// rot_pose( pose1_, pre_rotation_ );
		// rot_pose( pose2_, pre_rotation_ );

		clash_range_ = clash_range;
		inter_range_ = inter_range;
		grid_spacing_ = grid_spacing;
		min_contact_volume_ = min_contact_volume;
		aids_clash_1_ = aids_clash_1;
		aids_clash_2_ = aids_clash_2;
		aids_inter_1_ = aids_inter_1;
		aids_inter_2_ = aids_inter_2;
		filter_direction_ = filter_direction.normalized();
		filter_distance_ = filter_distance;

		numeric::xyzVector<double> bounds1;
		if( fixed ){
			// bounds1 = center_residues_rect( pose1_, aids_inter_1_, bounds1/*dummy*/, false );
			radius1 = center_cbs_welzl(pose1_,aids_inter_1_,center1_,true);
			bounds1 = numeric::xyzVector<double>(radius1,radius1,radius1);
			center1_ = 0;
			radius2 = center_cbs_welzl(pose2_,aids_inter_2_,center2_,true);
			center2_ = 0;
		} else {
			radius1 = center_cbs_welzl(pose1_,aids_inter_1_,center1_,false);
			bounds1 = numeric::xyzVector<double>(radius1,radius1,radius1);
			radius2 = center_cbs_welzl(pose2_,aids_inter_2_,center2_,false);
		}

		if(verbose_) std::cout << "FFTWDock: body welzl radii: " << radius1 << ", " << radius2 << std::endl;

		gsize_x_ = 1 + 2 * (core::Size)std::ceil( 2*std::max(bounds1[0],radius2) + inter_range_ ) / grid_spacing_ ;
		gsize_y_ = 1 + 2 * (core::Size)std::ceil( 2*std::max(bounds1[1],radius2) + inter_range_ ) / grid_spacing_ ;
		gsize_z_ = 1 + 2 * (core::Size)std::ceil( 2*std::max(bounds1[2],radius2) + inter_range_ ) / grid_spacing_ ;

		grid_clash_1_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		grid_clash_2_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		grid_inter_1_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		grid_inter_2_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		freq_clash_1_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		freq_clash_2_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		freq_inter_1_.dimension( gsize_x_, gsize_y_, gsize_z_ );
		freq_inter_2_.dimension( gsize_x_, gsize_y_, gsize_z_ );

		if( verbose_ ) std::cout << "FFTWDock: grid size: " << gsize_x_ <<" "<< gsize_y_ <<" "<< gsize_z_ << " ncells: " << gsize_x_*gsize_y_*gsize_z_ << std::endl;

		make_clash_grid( pose1_, aids_clash_1_, grid_clash_1_ );
		make_inter_grid( pose1_, aids_inter_1_, grid_inter_1_ );

		make_fftw_plans( gsize_z_, gsize_y_, gsize_x_ );
		plans_made_ = true;

		fft_forward( grid_clash_1_, freq_clash_1_ );
		fft_forward( grid_inter_1_, freq_inter_1_ );
	 }

	void get_hits_for_rotation( numeric::xyzMatrix<core::Real> const & rot1, numeric::xyzMatrix<core::Real> const & rot2, FFTHits & hits ){

		runtime_assert_msg( pre_rotation_ == numeric::xyzMatrix<core::Real>::identity(), "pre_rotation not implemented (correctly)" );
		numeric::xyzMatrix<core::Real> const r = pre_rotation_.transposed() * rot2 * pre_rotation_;
		rot_pose( pose2_, r );

		make_clash_grid( pose2_, aids_clash_2_, grid_clash_2_ );
		make_inter_grid( pose2_, aids_inter_2_, grid_inter_2_ );

				// dump_grid( grid_clash_1_, grid_spacing_, "grid_clash_1_.pdb" );
				// dump_grid( grid_clash_2_, grid_spacing_, "grid_clash_2_.pdb" );
				// dump_grid( grid_inter_1_, grid_spacing_, "grid_inter_1_.pdb" );
				// dump_grid( grid_inter_2_, grid_spacing_, "grid_inter_2_.pdb" );
				// pose1_.dump_pdb("pose1.pdb");
				// pose2_.dump_pdb("pose2.pdb");
				// utility_exit_with_message("dump grids");

	    fft_forward( grid_clash_2_, freq_clash_2_ );
	    fft_forward( grid_inter_2_, freq_inter_2_ );

	    // NOTE: re-use #2 grids for convolution!
	    for (int i=0; i<(int)freq_clash_2_.size(); ++i){
	        freq_clash_2_[i] = std::conj( freq_clash_1_[i] ) * freq_clash_2_[i]; // clash part
	        freq_inter_2_[i] = std::conj( freq_inter_1_[i] ) * freq_inter_2_[i]; // interaction part
	    }

		fft_reverse( freq_clash_2_, grid_clash_2_ );
		fft_reverse( freq_inter_2_, grid_inter_2_ );

		get_hits( grid_clash_2_, grid_inter_2_, rot1, rot2, hits );

	    // put back..
	    rot_pose( pose2_, r.transposed() );

	 }

	virtual void make_fftw_plans(int gsize_x, int gsize_y, int gsize_z) = 0;
	virtual void fft_forward( Complex3D & in, Complex3D & out ) = 0;
	virtual void fft_reverse( Complex3D & in, Complex3D & out ) = 0;
	virtual void get_hits(
		Complex3D const & clash_conv,
		Complex3D const & inter_conv,
		numeric::xyzMatrix<core::Real> const & rot1,
		numeric::xyzMatrix<core::Real> const & rot2,
		FFTHits & hits
	  ) = 0;

 protected:

	core::Real center_cbs_welzl(
		core::pose::Pose & pose,
		utility::vector1<core::id::AtomID> const & aids,
		numeric::xyzVector<core::Real> & delta,
		bool fixed = false
	 ){
		std::vector<numeric::xyzVector<core::Real> > cbs;
		for( auto i : aids ){
			cbs.push_back( pose.xyz(i) );
		}
		Sphere w = approx_bounding_sphere(cbs);
		if( !fixed ){
			for(core::Size ir = 1; ir <= pose.n_residue(); ++ir){
				for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia){
					core::id::AtomID aid(ia,ir);
					pose.set_xyz(aid,pose.xyz(aid)-w.center);
				}
			}
			delta = -w.center;
			return w.radius;
		} else {
			delta = 0;
			return w.radius + w.center.length();
		}
	 }
	numeric::xyzVector<core::Real> center_residues_rect(
		core::pose::Pose & pose,
		utility::vector1<core::id::AtomID> const & aids,
		numeric::xyzVector<core::Real> & delta,
		bool fixed = false
	 ){
		numeric::xyzVector<core::Real> lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);

		for( auto i : aids ){
			lb.min( pose.xyz(i) );
			ub.max( pose.xyz(i) );
		}
		if( !fixed ){
			numeric::xyzVector<core::Real> cen = (lb+ub)/2.0;
			for(core::Size ir = 1; ir <= pose.n_residue(); ++ir){
				for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia){
					core::id::AtomID aid(ia,ir);
					pose.set_xyz(aid,pose.xyz(aid)-cen);
				}
			}
			delta = -cen;
			return (ub-lb)/2.0;
		} else {
			return ub.max(-lb);
		}

	 }
	void make_clash_grid(
		core::pose::Pose const & pose,
		AIDs const & clash_aids,
		Complex3D & grid
	 ){
		grid = 0;

		numeric::xyzVector<core::Real> grid_offset = grid_spacing_ * ( numeric::xyzVector<core::Real>( grid.size3(), grid.size2(), grid.size1()  ) - 1.0 ) / 2.0;

		// set positive interactions
		for( auto aid : clash_aids ){
			numeric::xyzVector<core::Real> atom_xyz = pose.xyz(aid);
			core::Size lbx = (core::Size)std::max(            1.0           , (atom_xyz.x()+grid_offset.x() - clash_range_)/grid_spacing_ );
			core::Size lby = (core::Size)std::max(            1.0           , (atom_xyz.y()+grid_offset.y() - clash_range_)/grid_spacing_ );
			core::Size lbz = (core::Size)std::max(            1.0           , (atom_xyz.z()+grid_offset.z() - clash_range_)/grid_spacing_ );
			core::Size ubx = (core::Size)std::min( (core::Real)grid.size3() , (atom_xyz.x()+grid_offset.x() + clash_range_)/grid_spacing_ + 1.0 );
			core::Size uby = (core::Size)std::min( (core::Real)grid.size2() , (atom_xyz.y()+grid_offset.y() + clash_range_)/grid_spacing_ + 1.0 );
			core::Size ubz = (core::Size)std::min( (core::Real)grid.size1() , (atom_xyz.z()+grid_offset.z() + clash_range_)/grid_spacing_ + 1.0 );
			if( ubx > 2*grid_offset+1 ) continue;
			if( uby > 2*grid_offset+1 ) continue;
			if( ubz > 2*grid_offset+1 ) continue;
			for(core::Size x = lbx; x <= ubx; ++x){
			for(core::Size y = lby; y <= uby; ++y){
			for(core::Size z = lbz; z <= ubz; ++z){
				// std::cout << "make_grid x " << aid << lbx << " " << ubx << " " << x << std::endl;
				// std::cout << "make_grid y " << aid << lby << " " << uby << " " << y << std::endl;
				// std::cout << "make_grid z " << aid << lbz << " " << ubz << " " << z << std::endl;
				numeric::xyzVector<core::Real> cen = grid_spacing_*numeric::xyzVector<core::Real>(x,y,z) - grid_offset;
				if( cen.distance_squared( atom_xyz ) < clash_range_*clash_range_ ){
					grid(z,y,x) = 10000.0 * grid_spacing_ * grid_spacing_ * grid_spacing_ ;
				}
			}}}
		}

	 }
	void make_inter_grid(
		core::pose::Pose const & pose,
		AIDs const & inter_aids,
		Complex3D & grid
	 ){
		grid = 0;

		numeric::xyzVector<core::Real> grid_offset = grid_spacing_ * ( numeric::xyzVector<core::Real>( grid.size3(), grid.size2(), grid.size1()  ) - 1.0 ) / 2.0;

		// set positive interactions
		for( auto aid :inter_aids ){
			numeric::xyzVector<core::Real> atom_xyz = pose.xyz(aid);
			core::Size lbx = (core::Size)std::max(            1.0           , (atom_xyz.x()+grid_offset.x() - inter_range_)/grid_spacing_ );
			core::Size lby = (core::Size)std::max(            1.0           , (atom_xyz.y()+grid_offset.y() - inter_range_)/grid_spacing_ );
			core::Size lbz = (core::Size)std::max(            1.0           , (atom_xyz.z()+grid_offset.z() - inter_range_)/grid_spacing_ );
			core::Size ubx = (core::Size)std::min( (core::Real)grid.size3() , (atom_xyz.x()+grid_offset.x() + inter_range_)/grid_spacing_ + 1.0 );
			core::Size uby = (core::Size)std::min( (core::Real)grid.size2() , (atom_xyz.y()+grid_offset.y() + inter_range_)/grid_spacing_ + 1.0 );
			core::Size ubz = (core::Size)std::min( (core::Real)grid.size1() , (atom_xyz.z()+grid_offset.z() + inter_range_)/grid_spacing_ + 1.0 );
			for(core::Size x = lbx; x <= ubx; ++x){
			for(core::Size y = lby; y <= uby; ++y){
			for(core::Size z = lbz; z <= ubz; ++z){
				numeric::xyzVector<core::Real> cen = grid_spacing_*numeric::xyzVector<core::Real>(x,y,z) - grid_offset;
				if( cen.distance_squared( atom_xyz ) < inter_range_*inter_range_ ){
					grid(z,y,x) = sqrt( grid_spacing_ * grid_spacing_ * grid_spacing_ / (inter_range_ - clash_range_) );
				}
			}}}
		}

	 }
	inline int min_mod(int x,int y) {
		int r=x%y; if (r<-y/2) r+=y;if (r>=y/2) r-=y;
		return r;
	 }
	void dump_grid(
		Complex3D const & grid,
		double grid_spacing,
		std::string fname
	 ){
		numeric::xyzVector<core::Real> grid_offset = grid_spacing * ( numeric::xyzVector<core::Real>( grid.size3(), grid.size2(), grid.size1()  ) - 1.0 ) / 2.0;

		using namespace ObjexxFCL::format;
		utility::io::ozstream out(fname);
		core::Size i = 0;
		for(core::Size x = 1; x <= grid.size3(); ++x){
		for(core::Size y = 1; y <= grid.size2(); ++y){
		for(core::Size z = 1; z <= grid.size1(); ++z){
			if( grid(z,y,x).real() == 0.0) continue;
			std::string N = "  N ";
			if( fabs(grid(z,y,x).real() ) > 1.1 ) N = "  O ";
			++i;
			numeric::xyzVector<core::Real> viz = grid_spacing*numeric::xyzVector<core::Real>(x,y,z) - grid_offset;
			out<<"HETATM"<<I(5,i)<<' '<<N<<' ' <<	"VIZ"<<' '<<"B"<<I(4,i/10)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

		}}}

		out.close();
	 }
	void trans_pose( core::pose::Pose & pose, numeric::xyzVector<core::Real> const & trans ) {
		for(core::Size ir = 1; ir <= pose.n_residue(); ++ir) {
			for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
				core::id::AtomID const aid(core::id::AtomID(ia,ir));
				pose.set_xyz( aid, pose.xyz(aid) + trans );
			}
		}
	 }
	void rot_pose( core::pose::Pose & pose, numeric::xyzMatrix<core::Real> const & rot ) {
		for(core::Size ir = 1; ir <= pose.n_residue(); ++ir) {
			for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
				core::id::AtomID const aid(core::id::AtomID(ia,ir));
				pose.set_xyz( aid, rot * pose.xyz(aid) );
			}
		}
	 }

 };

struct FFTWDock3D : public FFTWDock {
	FFTWDock3D(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2
	) : FFTWDock(pose1,pose2) {}
	virtual ~FFTWDock3D(){}
	virtual void make_fftw_plans(int gsize_x, int gsize_y, int gsize_z){
		fftw_plan_forward_ = fftw_plan_dft_3d(
			gsize_x, gsize_y, gsize_z,
			(fftw_complex*)(&grid_clash_1_[0]),
	       	(fftw_complex*)(&freq_clash_1_[0]),
	       	FFTW_FORWARD,
	       	FFTW_MEASURE
	      );
		fftw_plan_reverse_ = fftw_plan_dft_3d(
			gsize_x, gsize_y, gsize_z,
			(fftw_complex*)(&grid_clash_1_[0]),
	       	(fftw_complex*)(&freq_clash_1_[0]),
	       	FFTW_BACKWARD,
	       	FFTW_MEASURE
	      );
	 }
	virtual void fft_forward( Complex3D & in, Complex3D & out ){
		runtime_assert( in.size3() == out.size3() );
		runtime_assert( in.size2() == out.size2() );
		runtime_assert( in.size1() == out.size1() );
	    fftw_execute_dft( fftw_plan_forward_, (fftw_complex*)(&in(1,1,1)), (fftw_complex*)(&out(1,1,1)) );
	 }
	virtual void fft_reverse( Complex3D & in, Complex3D & out ){
		runtime_assert( in.size3() == out.size3() );
		runtime_assert( in.size2() == out.size2() );
		runtime_assert( in.size1() == out.size1() );
	    fftw_execute_dft( fftw_plan_reverse_, (fftw_complex*)(&in(1,1,1)), (fftw_complex*)(&out(1,1,1)) );
	 }
	virtual void get_hits(
		Complex3D const & clash_conv,
		Complex3D const & inter_conv,
		numeric::xyzMatrix<core::Real> const & rot1,
		numeric::xyzMatrix<core::Real> const & rot2,
		FFTHits & hits
	){
		runtime_assert( clash_conv.size3() == inter_conv.size3() );
		runtime_assert( clash_conv.size2() == inter_conv.size2() );
		runtime_assert( clash_conv.size1() == inter_conv.size1() );
	    for(core::Size x = 1; x <= clash_conv.size3(); ++x){
	    for(core::Size y = 1; y <= clash_conv.size2(); ++y){
	    for(core::Size z = 1; z <= clash_conv.size1(); ++z){
	    	core::Real score = inter_conv(z,y,x).real() - clash_conv(z,y,x).real();
	    	score /= clash_conv.size(); // divide by N for fftw
	    	if( score < min_contact_volume_ ) continue;
			core::Real cx = grid_spacing_ * (core::Real)min_mod( x-1 , clash_conv.size3() ) ;
			core::Real cy = grid_spacing_ * (core::Real)min_mod( y-1 , clash_conv.size2() ) ;
			core::Real cz = grid_spacing_ * (core::Real)min_mod( z-1 , clash_conv.size1() ) ;
	    	numeric::xyzVector<core::Real> trans = -numeric::xyzVector<core::Real>(cx,cy,cz);
	    	if( fabs(trans.dot( filter_direction_ )) > filter_distance_ ) continue;

	    	// if( score > 439 ){
		    // 	core::pose::Pose tmp1(pose1_);
		    // 	core::pose::Pose tmp2(pose2_);
	    	// 	trans_pose(tmp1,center1_      );
	    	// 	trans_pose(tmp2,center2_+trans);
	    	// 	static int count = 0;
		    // 	tmp1.dump_pdb("FFTWDock_"+ObjexxFCL::format::I(4,++count)+"_1.pdb");
		    // 	tmp2.dump_pdb("FFTWDock_"+ObjexxFCL::format::I(4,  count)+"_2.pdb");
		    // 	std::cout << "FFTWDock DUMP " << score << " " << trans << std::endl;
		    // 	if(count > 100) std::exit(0);
	    	// }

	    	FFTHit hit;
	    	hit.score = score;
	    	hit.trans = pre_rotation_.transposed() * trans;
	    	hit.rot1 = rot1;
	    	hit.rot2 = rot2;
	    	hits.push_back(hit);
	    }}}

	 }

 };

struct FFTWDock2D_YZ : public FFTWDock {
	FFTWDock2D_YZ(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2
	) : FFTWDock(pose1,pose2) {}
	virtual ~FFTWDock2D_YZ(){}
	virtual void make_fftw_plans(int /*gsize_x*/, int gsize_y, int gsize_z){
		fftw_plan_forward_ = fftw_plan_dft_2d(
			gsize_y, gsize_z,
			(fftw_complex*)(&grid_clash_1_[0]),
	       	(fftw_complex*)(&freq_clash_1_[0]),
	       	FFTW_FORWARD,
	       	FFTW_MEASURE
	      );
		fftw_plan_reverse_ = fftw_plan_dft_2d(
			gsize_y, gsize_z,
			(fftw_complex*)(&grid_clash_1_[0]),
	       	(fftw_complex*)(&freq_clash_1_[0]),
	       	FFTW_BACKWARD,
	       	FFTW_MEASURE
	      );
	 }
	virtual void fft_forward( Complex3D & in, Complex3D & out ){
		runtime_assert( in.size3() == out.size3() );
		runtime_assert( in.size2() == out.size2() );
		runtime_assert( in.size1() == out.size1() );
		for(core::Size i = 1; i <= in.size3(); ++i){
		    fftw_execute_dft( fftw_plan_forward_, (fftw_complex*)(&in(1,1,i)), (fftw_complex*)(&out(1,1,i)) );
		}
	 }
	virtual void fft_reverse( Complex3D & in, Complex3D & out ){
		runtime_assert( in.size3() == out.size3() );
		runtime_assert( in.size2() == out.size2() );
		runtime_assert( in.size1() == out.size1() );
		for(core::Size i = 1; i <= in.size3(); ++i){
		    fftw_execute_dft( fftw_plan_reverse_, (fftw_complex*)(&in(1,1,i)), (fftw_complex*)(&out(1,1,i)) );
		}
	 }
	virtual void get_hits(
		Complex3D const & clash_conv,
		Complex3D const & inter_conv,
		numeric::xyzMatrix<core::Real> const & rot1,
		numeric::xyzMatrix<core::Real> const & rot2,
		FFTHits & hits
	){
		runtime_assert( clash_conv.size3() == inter_conv.size3() );
		runtime_assert( clash_conv.size2() == inter_conv.size2() );
		runtime_assert( clash_conv.size1() == inter_conv.size1() );
	    for(core::Size y = 1; y <= clash_conv.size2(); ++y){
	    for(core::Size z = 1; z <= clash_conv.size1(); ++z){
	    	core::Real score = 0;
		    for(core::Size x = 1; x <= clash_conv.size3(); ++x){
		    	score += inter_conv(z,y,x).real() - clash_conv(z,y,x).real();
			}
	    	score /= clash_conv.size() / clash_conv.size3(); // divide by N for fftw
	    	if( score < min_contact_volume_ ) continue;
			core::Real cx = 0; //grid_spacing_ * (core::Real)min_mod( x-1 , clash_conv.size3() ) ;
			core::Real cy = grid_spacing_ * (core::Real)min_mod( y-1 , clash_conv.size2() ) ;
			core::Real cz = grid_spacing_ * (core::Real)min_mod( z-1 , clash_conv.size1() ) ;
	    	numeric::xyzVector<core::Real> trans = -numeric::xyzVector<core::Real>(cx,cy,cz);
			if( fabs(trans.dot( filter_direction_ )) > filter_distance_ ) continue;

	    	// if( score > 700 ){
		    // 	core::pose::Pose tmp(scaffold_rot);
	    	// 	trans_pose(tmp,trans);
	    	// 	static int count = 0;
		    // 	tmp.dump_pdb("scaffold_dock_"+ObjexxFCL::format::I(4,irot)+"_"+ObjexxFCL::format::I(4,++count)+".pdb");
		    // 	std::cout << score << " " << trans << std::endl;
		    // 	if(count > 100) std::exit(0);
	    	// }

	    	FFTHit hit;
	    	hit.score = score;
	    	hit.trans = pre_rotation_.transposed() * trans;
	    	hit.rot1 = rot1;
	    	hit.rot2 = rot2;
	    	hits.push_back(hit);
	    }}

	 }

 };

#endif
