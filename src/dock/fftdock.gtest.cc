#include <gtest/gtest.h>

#include "dock/fftdock.hh"
#include <boost/multi_array.hpp>
#include "util/SimpleArray.hh"
#include <complex>

namespace scheme { namespace dock { namespace test {

using std::cout;
using std::endl;

TEST(fftdock,fftw){
	boost::multi_array<std::complex<double>,3> grid   (boost::extents[3][2][4]);
	boost::multi_array<std::complex<double>,3> grid_f (boost::extents[3][2][4]);	
	boost::multi_array<std::complex<double>,3> grid_ff(boost::extents[3][2][4]);	
	size_t N = grid.shape()[0] * grid.shape()[1] * grid.shape()[2];
	std::fill( grid   .origin(), grid   .origin() + grid   .size(), 0 );
	std::fill( grid_f .origin(), grid_f .origin() + grid_f .size(), 0 );
	std::fill( grid_ff.origin(), grid_ff.origin() + grid_ff.size(), 0 );		
	grid[1][0][2] = 1;

	fftw_plan p;
    p = fftw_plan_dft_3d( grid.shape()[0], grid.shape()[1], grid.shape()[2],
                         (fftw_complex*)grid.data(), (fftw_complex*)grid_f.data(), FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute(p);
    fftw_destroy_plan(p);
    p = fftw_plan_dft_3d( grid.shape()[0], grid.shape()[1], grid.shape()[2],
                         (fftw_complex*)grid_f.data(), (fftw_complex*)grid_ff.data(), FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute(p);
    fftw_destroy_plan(p);

    for(int i = 0; i < grid.shape()[0]; ++i)
    for(int j = 0; j < grid.shape()[1]; ++j)
    for(int k = 0; k < grid.shape()[2]; ++k)
    	ASSERT_NEAR( grid[i][j][k].real(), grid_ff[i][j][k].real()/N, 1e-9 );

 //    for(int i = 0; i < grid.shape()[0]; ++i){
 //    for(int j = 0; j < grid.shape()[1]; ++j){
 //    for(int k = 0; k < grid.shape()[2]; ++k){
 //    	cout << grid[i][j][k] << " ";
 //    } cout << endl;
 //    } cout << endl;
	// }
 //    for(int i = 0; i < grid.shape()[0]; ++i){
 //    for(int j = 0; j < grid.shape()[1]; ++j){
 //    for(int k = 0; k < grid.shape()[2]; ++k){
 //    	cout << grid_f[i][j][k] << " ";
 //    } cout << endl;
 //    } cout << endl;
	// }
 //    for(int i = 0; i < grid.shape()[0]; ++i){
 //    for(int j = 0; j < grid.shape()[1]; ++j){
 //    for(int k = 0; k < grid.shape()[2]; ++k){
 //    	cout << grid_ff[i][j][k] << " ";
 //    } cout << endl;
 //    } cout << endl;
	// }
}

}}}
