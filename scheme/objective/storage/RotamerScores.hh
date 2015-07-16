#ifndef INCLUDED_objective_storage_RotamerScores_HH
#define INCLUDED_objective_storage_RotamerScores_HH

#include "scheme/util/SimpleArray.hh"

#include <boost/lexical_cast.hpp>

namespace scheme { namespace objective { namespace storage {

template< class _Data = uint16_t, int _RotamerBits = 9, int _Divisor = -13 >
struct RotamerScore {
	typedef RotamerScore< _Data, _RotamerBits, _Divisor > THIS;
	typedef _Data Data;
	static const int RotamerBits = _RotamerBits;
	static const int ScoreBits = sizeof(_Data)*8 - _RotamerBits;	
	static const Data one = 1;
	static const Data RotamerMask = (( one << RotamerBits ) - one );

	Data data_;
	RotamerScore(): data_( RotamerMask ) {}
	RotamerScore( Data data ): data_(data) {}
	RotamerScore( Data rot, float score ){
		Data sdat = score * _Divisor;
		assert( sdat < ( one<<ScoreBits ) );
		assert( rot < (one<<RotamerBits) );
		data_ = rot | (sdat<<RotamerBits);
	}
	float score() const { return data2float( get_score_data() ); }
	Data  rotamer() const { return data_ & RotamerMask; }
	void  set_score( float score ){ set_score_data( float2data( score ) ); }
	void  set_rotamer( Data rot ){ assert(rot < one<<RotamerBits); data_ = (data_& (~RotamerMask)) | rot; }

	void  set_score_data( Data sd ){ assert(sd < (one<<ScoreBits)); data_ = rotamer() | ( sd << RotamerBits ); }
	Data  get_score_data() const { return data_ >> RotamerBits; }
	static float divisor() { return _Divisor; }

	static float data2float( Data data ){ return float(data)/_Divisor; }
	static Data  float2data( float f ){ return Data( f*_Divisor ); }

	bool operator < ( THIS const & other ) const { return data_ > other.data_; } // reverse so low score is low
	bool operator== ( THIS const & other ) const { return data_ == other.data_; }
	bool operator== ( Data const & other ) const { return data_ == other; }	

	bool empty() const { return data_ == RotamerMask; }
};

template< int _N, class _Data = uint16_t, int _RotamerBits=9, int _Divisor=-13 >
struct RotamerScores {
	BOOST_STATIC_ASSERT( _N > 0   );
	BOOST_STATIC_ASSERT( _N < 256 ); // arbitrary

	typedef _Data Data;
	typedef RotamerScores< _N, _Data, _RotamerBits, _Divisor > THIS;
	typedef RotamerScore<      _Data, _RotamerBits, _Divisor > RotScore;
	
	static int const N = _N;
	util::SimpleArray<N,RotScore> rotscores_;
	
	RotamerScores(){
		rotscores_.fill( RotScore::RotamerMask );
	}
	
	void add_rotamer( Data rot, float score ){
		add_rotamer( RotScore(rot,score) );
	}
	void add_rotamer( RotScore to_insert ){
		Data irot = to_insert.rotamer();
		int insert_pos = 0;
		RotScore worst( std::numeric_limits<Data>::max() );
		for( int i = 0; i < N; ++i ){
			// std::cout << " iter " << i << " " << rotscores_[i].score() << " " << rotscores_[i].rotamer() << " " 
			          // << "cur " << rotscores_[i].data_ << " low " << worst.data_ << std::endl;
			// if rot already stored, this is the position we check
			if( rotscores_[i].rotamer() == irot ){
				insert_pos = i;
				// std::cout << "rotamer equal at " << i << std::endl;
				break;
			}
			// else we take the worst position
			if( worst < rotscores_[i] ){
				// std::cout << "worst is " << i << std::endl;
				worst = rotscores_[i];
				insert_pos = i;
			}
		}
		// std::cout << "insert_pos " << insert_pos << std::endl;
		// now insert if new val is better than worst stored val
		if( to_insert < rotscores_[insert_pos] ){
			rotscores_[insert_pos] = to_insert;
		}
	}

	void merge( THIS const & other ){
		int size = other.size();
		for( int i = 0; i < size; ++i ){
			add_rotamer( other.rotscores_[i] );
		}
	}

	float score_of_rotamer( int irot ) const {
		for( int i = 0; i < N; ++i ){
			if( rotscores_[i].rotamer() == irot ){
				return rotscores_[i].score();
			}
		}
		return 0.0f;
	}

	float score( int irot ) const { assert(irot<N); return rotscores_[irot].score(); }
	float rotamer( int irot ) const { assert(irot<N); return rotscores_[irot].rotamer(); }	

	static int maxsize(){ return _N; }

	int size() const { int i; for(i=0;i<_N;++i) if( rotscores_[i].empty() ) break; return i; }	

	void sort_rotamers(){
		std::sort( rotscores_.begin(), rotscores_.end() );
	}

	bool is_sorted() const {
		for( int i = 1; i < _N; ++i )
			if( rotscores_[i] < rotscores_[i-1] )
				return false;		
		return true;
	}

	static std::string name() {
		static std::string const name = "RotamerScores< "
			 + boost::lexical_cast<std::string>(_N)
			 + ", "
			 + boost::lexical_cast<std::string>(sizeof(_Data))	 
			 + ", "
			 + boost::lexical_cast<std::string>(_RotamerBits)
			 + ", "
			 + boost::lexical_cast<std::string>(_Divisor)
			 +" >";
		return name;
	}

};

template< int N, class V, int B, int D >
std::ostream & operator << ( std::ostream & out, RotamerScores<N,V,B,D> const & val ){
	out << "RotamerScores( ";
	for(int i = 0; i < val.size(); ++i){
		out << val.score(i) << "," << val.rotamer(i) << "  ";
	}
	out << ")";
	return out;
}

}}}

#endif
