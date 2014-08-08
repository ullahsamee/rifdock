#ifndef INCLUDED_io_cache_HH
#define INCLUDED_io_cache_HH

#include <iostream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


namespace scheme { namespace io {

	template<class T>
	bool read_cache(std::string const & location, T & t){
		using namespace boost::iostreams;
		if( location.size()==0 || !boost::filesystem::exists(location) ) return false;
		// std::cout << "read_cache " << location << std::endl;
		std::ifstream file( location.c_str() , std::ios_base::in|std::ios_base::binary );
		filtering_streambuf<input> in;
		in.push(gzip_decompressor());
		in.push(file);
		std::istream stdin(&in);
        boost::archive::binary_iarchive archive(stdin);
        archive >> t;
		return true;
	}

	template<class T>
	void write_cache(std::string const & location, T & t){
		if(location.size()==0) return;
		// std::cout << "write_cache " << location << std::endl;
		using namespace boost::iostreams;
		std::ofstream file( location.c_str() , std::ios_base::out|std::ios_base::binary );
		if( !file.good() ) throw 1;
		filtering_streambuf<output> out;
		out.push(gzip_compressor());
		out.push(file);
		std::ostream stdout(&out);
        boost::archive::binary_oarchive archive(stdout);
        archive << t;
	}

}}

#endif
