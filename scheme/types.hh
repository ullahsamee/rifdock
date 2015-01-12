#ifdef CXX11
	#include <memory>
#else
	#include <boost/shared_ptr.hpp>
	#include <boost/make_shared.hpp>
	#include <boost/enable_shared_from_this.hpp>
#endif

namespace scheme {
	
	#ifdef CXX11
		using std::shared_ptr;
		using std::make_shared;
		using std::enable_shared_from_this;
	#else
		using boost::shared_ptr;
		using boost::make_shared;
		using boost::enable_shared_from_this;
	#endif


}
