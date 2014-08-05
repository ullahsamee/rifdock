#ifndef INCLUDED_kinematics_Scene_io_HH
#define INCLUDED_kinematics_Scene_io_HH

#include "kinematics/Scene.hh"

#include "util/meta/print_type.hh"

#include <boost/foreach.hpp>

namespace scheme {
namespace kinematics {

namespace m = boost::mpl;
namespace f = boost::fusion;
using std::cout;
using std::endl;


	template< class C, class P, class I >
	std::ostream & operator<<(std::ostream & out, impl::BodyTplt<C,P,I> const & b){
		return out << "Body( " << b.position() << ", " << &b.conformation() << " )";
	}

	template<class C,class P>
	std::ostream & operator<<(std::ostream & out, Scene<C,P> const & scene){ 
		using std::endl;
		typedef Scene<C,P> Scene;
		out << "Scene" << endl;
		out << "  Nbodies: Asym: " << scene.nbodies_asym() << " Total: " << scene.nbodies() << endl;
		BOOST_FOREACH(typename Scene::Body const & b, scene.__bodies_unsafe__()){ out << "    " << b << endl; }
		out << "  Types:" << endl;
		out << "    Actors:" << endl; util::meta::print_type<typename Scene::Actors>(out,"        ");
		out << "    Conformation:" << endl; util::meta::print_type<typename Scene::Conformation>(out,"        ");	
		return out;
	}

}
}

#endif
