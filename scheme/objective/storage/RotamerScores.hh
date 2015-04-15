#ifndef INCLUDED_objective_storage_RotamerScores_HH
#define INCLUDED_objective_storage_RotamerScores_HH

#include "scheme/util/SimpleArray.hh"

#include <boost/lexical_cast.hpp>

namespace scheme { namespace objective { namespace storage {



template< int _N, class _Value = int16_t, int _Divisor=1024 >
struct RotamerScores {
	BOOST_STATIC_ASSERT( _N > 0   );
	BOOST_STATIC_ASSERT( _N < 256 ); // arbitrary

	typedef RotamerScores<_N,_Value,_Divisor> THIS;
	typedef _Value Value;
	
	static int const N = _N;
	static int const Divisor = _Divisor;
	static float const MUL =      _Divisor;
	static float const DIV = 1.0f/_Divisor;	
	
	util::SimpleArray<N,Value> rotamer_, scores_;
	
	RotamerScores(){
		rotamer_.fill( std::numeric_limits<Value>::max() );
		scores_ .fill( std::numeric_limits<Value>::max() );
	}
	
	void add_rotamer( Value rot, float sc ){
		Value score = float2val(sc);
		int   iworst = -1;
		Value vworst = std::numeric_limits<Value>::min();
		for( int i = 0; i < N; ++i){
			if( rotamer_.at(i) == rot ){ // if rot already stored, use it
				vworst = scores_.at(i);
				iworst = i;
				break;				
			}
			if( vworst <= scores_.at(i) ){ // else pick worst score
				vworst = scores_.at(i);
				iworst = i;
			}
		}
		// std::cout << "iworst " << iworst << " " << vworst << " " << score <<  std::endl;
		if( iworst < 0 || vworst <= score ) return;
		scores_ [iworst] = score;
		rotamer_[iworst] = rot;
	}

	float score_of_rot( int16_t rot ) const {
		Value sc = 0;
		for( int i = 0; i < N; ++i ){
			sc = ( rot == rotamer_.at(i) ) ? scores_.at(i) : sc;
			// std::cout << "score_of_rot " << i << " " << sc << " " << val2float(sc) << " DIV " << DIV << " " << float(sc) << " " << float(sc)*DIV << std::endl;
		}
		return val2float(sc);
	}

	float score( int i ) const { return val2float(scores_.at(i)); }
	int rotamer( int i ) const { return rotamer_.at(i); }
	
	static Value float2val( float f ) { return Value(f *MUL); }
	static float val2float( Value v ) { return float(v)*DIV ; }

	static int size(){ return _N; }

	void sort_rotamers(){
		for(int i = 0; i < N; ++i){
			for( int j = i+1; j < N; ++j){
				if( scores_[j] < scores_[i] ){
					std::swap( scores_ [j], scores_ [i] );
					std::swap( rotamer_[j], rotamer_[i] );					
				}
			}
		}
	}

	static std::string name() {
		static std::string const name = "RotamerScores< "
			 + boost::lexical_cast<std::string>(_N)
			 + ", "
			 + boost::lexical_cast<std::string>(sizeof(_Value))	 
			 + ", "
			 + boost::lexical_cast<std::string>(_Divisor)
			 +" >";
		return name;
	}

};



}}}

#endif
