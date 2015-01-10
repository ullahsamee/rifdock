#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
// #include <memory>

namespace scheme {
	
	#ifndef CXX11
		using boost::shared_ptr;
		using boost::make_shared;
		using boost::enable_shared_from_this;
	#else
		using std::shared_ptr;
		using std::make_shared;
		using std::enable_shared_from_this;
	#endif


}
