#ifndef HOUTILS_H
#define HOUTILS_H

#include "hoNullModel.h"
#include "../refinementPhase/MotifContainer.h"

using std::cout;
using std::endl;

#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif


double addLogspace( double x, double y );

std::vector<float> alphaPseudoCountsFactor( std::vector<float> alpha, float q );

char* baseFileName( char* name );

void calculateWeights( std::vector<float> &x, float x_bg, bool ranks=false );


template<typename T>
void printInterpolatedMarkovBackgroundModel( const T a, const int fields,
		                                     const int* breaks );

template<typename T>
void printInterpolatedMarkovModel( const T a, const int fr, const int to,
		                           const int fields, const int* breaks );
void printInterpolatedMarkovModel( Motif& model, bool parameters=false );
void printInterpolatedMarkovModels( MotifContainer& models,
		                            bool parameters=false );

float quantile( std::vector<float> x, float prob, int type=7 );
std::vector<float> rank( std::vector<float> x, bool decreasing=false );

void saveInterpolatedMarkovModel( Motif& model, const char* baseFileName );
void saveInterpolatedMarkovModels( MotifContainer& models,
		                           const char* baseFileName );

/* log( exp( x ) + exp( y ) )
 *   = x + log( 1 + exp( y-x ) )
 * x and y in log space */
inline double addLogspace( double x, double y ){

    assert( !( x != x ) && "addLogspace: x is NaN." );
    assert( !( y != y ) && "addLogspace: y is NaN." );

    if( x <= log( DBL_MIN ) ) {
        return y;
    } else if( y <= log( DBL_MIN ) ) {
        return x;
    } else if( x < y ){
        return( y + log( 1.0 + exp( x-y ) ) );
    } else{
        return( x + log( 1.0 + exp( y-x ) ) );
    }
}

inline char* baseFileName( char* name ){

	int pos = 0, startPos = 0, lastPos = 0;

	// '.' index
	while( name[++pos] != '\0' ){
		if( name[pos] == '.' ){
			lastPos = pos-1;
		}
	}

	// '/' index
	while( --pos != 0 && name[pos] != '/' ){
		;
	}
	if( pos != 0 ){
		startPos = pos+1;
	}

	char* basename = ( char* )malloc( ( lastPos-startPos+2 )*sizeof( char ) );

	for( pos=startPos; pos <= lastPos; pos++ ){
		basename[ pos-startPos ] = name[pos];
	}
	basename[ pos-startPos ] = '\0';

	return basename;
}

inline void calculateWeights( std::vector<float> &x, float x_bg, bool ranks ){

	if( ranks ){
		/*
		 * rank-based weights
		 */

		float N = static_cast<float>( x.size() );
		std::vector<float> r = rank( x, true );

		float N1 = N + 1;
		float N1xbg = N1 - x_bg;
		float denominator = N - N1xbg;

		for( unsigned i=0; i < r.size(); i++ ){
			if( r.at( i ) < x_bg ){
				x.at( i ) = ( ( N1 - r.at( i ) ) - N1xbg ) / denominator;
			} else{
				x.at( i ) = 0.0f;
			}
		}
	} else{
		/*
		 * intensity-based weights
		 */

		float denominator = *std::max_element( x.begin(), x.end() ) - x_bg;
		for( unsigned i=0; i < x.size(); i++ ){
			if( x.at( i ) > x_bg ){
				x.at( i ) = ( x.at( i ) - x_bg ) / denominator;
			} else{
				x.at( i ) = 0.0f;
			}
		}
	}
}

template<typename T>
inline void printInterpolatedMarkovBackgroundModel( const T a, const int fields,
		                                            const int* breaks ){

	int i, k;

	for( k=1, i=1; i <= fields; ++i ){
		if( i==breaks[k] ){

			cout << endl; ++k;
		}
		cout << a[i] << " ";
	}
	cout << endl;
}

template<typename T>
inline void printInterpolatedMarkovModel( const T a, const int from,
		                                  const int to, const int fields,
		                                  const int* breaks ){

	int i, k, pos;

	for( pos=from; pos <= to; ++pos ){
		for( k=1, i=1; i <= fields; ++i ){
			if( i==breaks[k] ){

				cout << endl; ++k;
			}
			cout << a[pos][i] << " ";
		}
		cout << endl;
		if( !( pos+1 > to ) ) cout << endl;
	}
}

inline void printInterpolatedMarkovModel( Motif& model, bool parameters ){

	int from = model.getFirstMotifColumn();
	int to = model.getLastMotifColumn();

	int order = model.getOrder();
	int* offsets = model.getOffsets();

	if( parameters ){
		printf( " __________\n"
				"|          |\n"
				"| SETTINGS |\n"
				"|__________|\n\n" );

		cout << "order" << "\t\t\t" << order << endl;

		cout << "specificity factor q" << "\t" << model.getSpecificityFactor()
			 << endl;

		cout << "pseudo-counts factor a" << "\t";
		for( int i=0; i < order+1; i++ ){
			cout << model.getAlpha()[i] << " ";
		}
		cout << endl;

	}

	printf( " ________\n"
			"|        |\n"
			"| COUNTS |\n"
			"|________|\n\n" );
	printInterpolatedMarkovModel( model.getCounts(), from, to,
			                      offsets[order+1]-1, offsets );
	if( order > 0 ){
		printf( " _________\n"
				"|         |\n"
				"| COUNTSX |\n"
				"|_________|\n\n" );
		printInterpolatedMarkovModel( model.getCountsx(), from, to-1,
				                      offsets[order]-1, offsets );
	}
	printf( " ___________________________\n"
			"|                           |\n"
			"| CONDITIONAL PROBABILITIES |\n"
			"|___________________________|\n\n" );
	printInterpolatedMarkovModel( model.getConds(), from, to,
			                      offsets[order+1]-1, offsets );
	printf( " _____________________\n"
			"|                     |\n"
			"| TOTAL PROBABILITIES |\n"
			"|_____________________|\n\n" );
	printInterpolatedMarkovModel( model.getPWM(), from, to, offsets[order+1]-1,
			                      offsets );
}

inline void printInterpolatedMarkovModels( MotifContainer& models,
		                                   bool parameters ){

	list<Motif*>::const_iterator iter;

	int no = 1;
	for( iter=models.getMotifs().begin(); iter != models.getMotifs().end();
		 ++iter, ++no ){
		printf( " _______\n"
				"|       |\n"
				"| MODEL |\n"
				"|_______|\n"
				"\n"
				"no. %d\n\n", no );
		printInterpolatedMarkovModel( **iter, parameters );
	}
}

inline float quantile( std::vector<float> x, float prob, int type ){

	/*
	 * Estimating quantiles
	 *
	 * http://en.wikipedia.org/wiki/Quantile
	 */

	if( prob < 0.0f || prob > 1.0f ){
		fprintf( stderr, "Calculate quantile using prob = [0,1]\n" );
		exit(1);
	}

	sort( x.begin(), x.end() );
	int N = static_cast<int>( x.size() );

    switch ( type ){
	case 3 :
		if( prob <= ( 1.0f / 2.0f ) / static_cast<float>( N ) ){
			return( x.front() );
		} else{
			float h = static_cast<float>( ( N ) )*prob;
			return( x.at( static_cast<int>( roundf( h ) )-1 ) );
		}
	case 7 :
		if( prob == 1.0f ){
			return( x.back() );
		} else{
			float h = static_cast<float>( ( N-1 ) )*prob + 1.0f;
			int h_ = static_cast<int>( floorf( h ) );
			return( x.at( h_-1 ) + ( h - static_cast<float>( h_ ) )*( x.at( h_ )
					- x.at( h_-1 ) ) );
		}
	case 8 :
		if( prob < ( 2.0f / 3.0f ) / ( static_cast<float>( N ) + ( 1.0f / 3.0f )
			) ){
			return( x.front() );
		} else if( prob >= ( static_cast<float>( N ) - ( 1.0f / 3.0f ) ) / (
				   static_cast<float>( N ) + ( 1.0f / 3.0f ) ) ){
			return( x.back() );
		} else{
			float h = ( static_cast<float>( N ) + ( 1.0f / 3.0f ) )*prob + (
					    1.0f / 3.0f );
			int h_ = static_cast<int>( floorf( h ) );
			return( x.at( h_-1 ) + ( h - static_cast<float>( h_ ) )*( x.at( h_ )
					- x.at( h_-1 ) ) );
		}
	default :
		fprintf( stderr, "Calculate quantile using type = {3, 7,8}\n" );
		exit( EXIT_FAILURE );
    }
}

inline std::vector<float> rank( std::vector<float> x, bool decreasing ){

  std::vector<int> idx( x.size() );
  for( int i=0; i < static_cast<int>( idx.size() ); i++ ){
	   idx[i] = i;
  }

  if( decreasing ){
	  sort( idx.begin(), idx.end(),
			[&x]( int i1, int i2 ){ return x[i1] > x[i2]; } );
  } else{
	  sort( idx.begin(), idx.end(),
			[&x]( int i1, int i2 ){ return x[i1] < x[i2]; } );
  }

  std::vector<float> ranks( x.size() );
  for( int i=0; i < static_cast<int>( idx.size() ); i++ ){
	   ranks[idx[i]] = static_cast<float>( i+1 );
  }

  return( ranks );
}

inline void saveInterpolatedMarkovModel( Motif& model, const char* baseFileName
		                                 ){
	/*
	 * Save interpolated Markov model to three flat files
	 * (1) baseFileName.freqs (background frequencies)
	 * (2) baseFileName.probs (probabilities)
	 * (3) baseFileName.conds (conditional probabilities)
	 */

	std::stringstream str;
	FILE *f_frequencies;
	FILE *f_conditionals;
	FILE *f_probabilities;

	// save background frequencies

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".freqs";
	f_frequencies = fopen( str.str().c_str(), "w" );

	for( int i=1; i <= nAlpha( Global::A ); ++i ){
		fprintf( f_frequencies, "%.8e ", Global::posBg[i] );
	}
	fprintf( f_frequencies, "\n" );
	fclose( f_frequencies );

	// save probabilities/conditional probabilities

	double** conds = model.getConds();
	double** probs = model.getPWM();
	int* offsets = model.getOffsets();
	int fields = offsets[model.getOrder()+1]-1;

	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".conds";
	f_conditionals = fopen( str.str().c_str(), "w" );

	str.str( "" );
	str<< Global::outputDirectory << "/" << baseFileName << ".probs";
	f_probabilities = fopen( str.str().c_str(), "w" );

	for( int pos=model.getFirstMotifColumn(); pos <= model.getLastMotifColumn();
		 ++pos ){
		for( int k=1, i=1; i <= fields; ++i ){
			if( i==offsets[k] ){
				fprintf( f_conditionals, "\n" );
				fprintf( f_probabilities, "\n" );
				++k;
			}
			fprintf( f_conditionals, "%.8e ", conds[pos][i] );
			fprintf( f_probabilities, "%.8e ", probs[pos][i] );
		}
		fprintf( f_conditionals, "\n\n" );
		fprintf( f_probabilities, "\n\n" );
	}
	fclose( f_conditionals );
	fclose( f_probabilities );
}

inline void saveInterpolatedMarkovModels( MotifContainer& models,
		                                  const char* baseFileName ){

	/*
	 * Save single interpolated Markov models to three flat files
	 * (1) baseFileName.freqs (background frequencies)
	 * (2) baseFileName.probs (probabilities)
	 * (3) baseFileName.conds (conditional probabilities)
	 */

	std::stringstream str;
	list<Motif*>::const_iterator iter;

	if( models.getMotifs().size() > 1 ){
		int no = 1; // number models
		for( iter=models.getMotifs().begin(); iter != models.getMotifs().end();
			 ++iter, ++no ){
			str.str( "" );
			str << baseFileName << "-" << no;
			saveInterpolatedMarkovModel( **iter, str.str().c_str() );
		}
	} else{
		saveInterpolatedMarkovModel( *models.getMotifs().front(), baseFileName
				                     );
	}
}





//
/** The digamma function in long double precision.
* @param x the real value of the argument
* @return the value of the digamma (psi) function at that point
* @author Richard J. Mathar
* @since 2005-11-24
*
* this piece of code is needed for the gradient calculation
* within the conjugate gradient for the alpha learning
*/
inline long double digamma(long double x)
{
	/* force into the interval 1..3 */
	if( x < 0.0L )
		return digamma(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
	else if( x < 1.0L )
		return digamma(1.0L+x)-1.0L/x ;
	else if ( x == 1.0L)
		return -M_GAMMAl ;
	else if ( x == 2.0L)
		return 1.0L-M_GAMMAl ;
	else if ( x == 3.0L)
		return 1.5L-M_GAMMAl ;
	else if ( x > 3.0L)
		/* duplication formula */
		return 0.5L*(digamma(x/2.0L)+digamma((x+1.0L)/2.0L))+M_LN2l ;
	else
	{
		/* Just for your information, the following lines contain
		* the Maple source code to re-generate the table that is
		* eventually becoming the Kncoe[] array below
		* interface(prettyprint=0) :
		* Digits := 63 :
		* r := 0 :
		*
		* for l from 1 to 60 do
		* 	d := binomial(-1/2,l) :
		* 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
		* 	evalf(r) ;
		* 	print(%,evalf(1+Psi(1)-r)) ;
		*o d :
		*
		* for N from 1 to 28 do
		* 	r := 0 :
		* 	n := N-1 :
		*
 		*	for l from iquo(n+3,2) to 70 do
		*		d := 0 :
 		*		for s from 0 to n+1 do
 		*		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
 		*		od :
 		*		if 2*l-n > 1 then
 		*		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
 		*		fi :
 		*	od :
 		*	print(evalf((-1)^n*2*r)) ;
 		*od :
 		*quit :
		*/
		static long double Kncoe[] = { .30459198558715155634315638246624251L,
		.72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
		.27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
		.17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
		.11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
		.83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
		.59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
		.42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
		.304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
		.21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
		.15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
		.11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
		.80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
		.58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
		.41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;

		register long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
		register long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
		register long double resul = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0L ;

		for(unsigned int n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
		{
			const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			resul += Kncoe[n]*Tn1 ;
			Tn_1 = Tn ;
			Tn = Tn1 ;
		}
		return resul ;
	}
}


inline void calcQandGrad(Motif &motif, unsigned char* kmer, int order, double alpha){
	int y,p, modelPos;
	list<int>::iterator j;
	double** conds = motif.getConds();
	double** counts = motif.getCounts();

	for(j = motif.getMotifColumns().begin(); j != motif.getMotifColumns().end(); j++){

		double Qfunc_p2_pos = 0.0;
		double Qfunc_p3_pos = 0.0;
		double Qfunc_p4_pos = 0.0;

		double Qfunc_grad_pos = 0.0;
		double Grad_p2_pos    = 0.0;
		double Grad_p3_pos    = 0.0;

		double zeroOrder_pos   = 0.0;
		double pseudoOrder_pos = 0.0;

		modelPos = *j - *motif.getMotifColumns().begin(); //0,1,...
		if( modelPos < order ){ // margins
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// longer kmer indices
				y = motif.sub2ind( kmer, 0, order );
				// pseudo-counts kmer indices
				p = motif.sub2ind(kmer+(order-modelPos), 0, modelPos);

				if(order == 0){
					// calculate contribution to Qfunc
					Qfunc_p2_pos += counts[*j][y] * log(conds[*j][y]);

//					printf("conds[%d][%d] = %f\n",*j,y,conds[*j][y]);

					Qfunc_p3_pos -= lgamma(alpha * Global::posBg[k] +1);
					Qfunc_p4_pos += Global::posBg[k] * log(conds[*j][y]);

					// calculate contribution to Gfunc_grad
					Qfunc_grad_pos += Global::posBg[k] * (double(digamma(alpha * Global::posBg[k] + 1))- log(conds[*j][y]));
					Grad_p2_pos    += Global::posBg[k] * double(digamma(alpha*Global::posBg[k] +1));
					Grad_p3_pos    += Global::posBg[k] * log(conds[*j][y]);

					zeroOrder_pos   += counts[*j][y] * log(conds[*j][y]);
					pseudoOrder_pos += counts[*j][y] * log(Global::posBg[k]);
				}else{
					// calculate contribution to Qfunc
					Qfunc_p2_pos += counts[*j][y] * log(conds[*j][y]);

//					printf("conds[%d][%d] = %f\n",*j,y,conds[*j][y]);

					Qfunc_p3_pos -= lgamma(alpha * conds[*j][p] +1);
					Qfunc_p4_pos += conds[*j][p] * log(conds[*j][y]);

					// calculate contribution to Gfunc_grad
					Qfunc_grad_pos += conds[*j][p] * (double(digamma(alpha * conds[*j][p] + 1))- log(conds[*j][y]));
					Grad_p2_pos    += conds[*j][p] * double(digamma(alpha*conds[*j][p] +1));
					Grad_p3_pos    += conds[*j][p] * log(conds[*j][y]);
				}

			} // end for loop over k
		}//end if
		else{ // remaining model positions
			for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
				// last kmer base
				kmer[order] = k;
				// longer kmer indices
				y = motif.sub2ind( kmer, 0,order );
				// pseudo-counts kmer indices
				p = motif.sub2ind( kmer+1,0,order-1);

				if(order == 0){
					// calculate contribution to Qfunc
					Qfunc_p2_pos += counts[*j][y] * log(conds[*j][y]);

//					printf("conds[%d][%d] = %f\n",*j,y,conds[*j][y]);

					Qfunc_p3_pos -= lgamma(alpha * Global::posBg[k] +1);
					Qfunc_p4_pos += Global::posBg[k] * log(conds[*j][y]);

					// calculate contribution to Gfunc_grad
					Qfunc_grad_pos += Global::posBg[k] * (double(digamma(alpha * Global::posBg[k] + 1))- log(conds[*j][y]));
					Grad_p2_pos    += Global::posBg[k] * double(digamma(alpha*Global::posBg[k] +1));
					Grad_p3_pos    += Global::posBg[k] * log(conds[*j][y]);

					zeroOrder_pos   += counts[*j][y] * log(conds[*j][y]);
					pseudoOrder_pos += counts[*j][y] * log(Global::posBg[k]);

				}else{
					// calculate contribution to Qfunc
					Qfunc_p2_pos += counts[*j][y] * log(conds[*j][y]);

//					printf("conds[%d][%d] = %f\n",*j,y,conds[*j][y]);

					Qfunc_p3_pos -= lgamma(alpha * conds[*j][p] +1);
					Qfunc_p4_pos += conds[*j][p] * log(conds[*j][y]);

					// calculate contribution to Gfunc_grad
					Qfunc_grad_pos += conds[*j][p] * (double(digamma(alpha * conds[*j][p] + 1))- log(conds[*j][y]));
					Grad_p2_pos    += conds[*j][p] * double(digamma(alpha*conds[*j][p] +1));
					Grad_p3_pos    += conds[*j][p] * log(conds[*j][y]);
				}

			}//end for loop over k
		} //end else
		motif.setQfunc_p1(lgamma(alpha + 4));
		motif.setQfunc_p2(Qfunc_p2_pos);
		motif.setQfunc_p3(Qfunc_p3_pos);
		motif.setQfunc_p4(Qfunc_p4_pos *alpha);

		motif.setQfunc_grad(Qfunc_grad_pos);
		motif.setGrad_p2(Grad_p2_pos);
		motif.setGrad_p3(Grad_p3_pos);

		motif.setzeroOrder(zeroOrder_pos);
		motif.setpseudoOrder(pseudoOrder_pos);
	} // end of loop j=0 -> W-1
}


inline void calcQandGrad(Motif &motif, unsigned char* kmer, int order, double alpha, int pos){
	if(order == 0){
    	calcQandGrad(motif,kmer,order, alpha);
	}
	else{
	  for( unsigned char k=1; k <= nAlpha( Global::A ); k++ ){
	    kmer[pos] = k; // next kmer base
	    if( pos == ( order-1 )){
	    	// calculate score
	    	calcQandGrad(motif,kmer,order, alpha);
	    } else{
		// add another base
		calcQandGrad(motif,kmer,order, alpha,pos+1);
	    }
	  }
	}
}



#endif /* HOUTILS_H */
