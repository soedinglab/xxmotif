/*
 * UniversalKmer.h
 *
 *  Created on: Oct 23, 2009
 *      Author: eckhart
 */

#ifndef UNIVERSALKMER_H_
#define UNIVERSALKMER_H_

#include "AbstractKmer.h"

#include <stdint.h>
#include <string>
#include <sstream>
#include <list>
#include <vector>

class UniversalKmer: public AbstractKmer {

	friend class UniversalKmerTest;
	//FRIEND_TEST(UniversalKmerTest, MutateAndLengthTest);

	private:
		  void operator=(const UniversalKmer&);

		class rep_type {
			public:
				char chr;
				int gap;
				rep_type(const char c, const int g) :
					chr(c), gap(g) {
				}
		};

		typedef std::vector<rep_type> list_type;
		list_type cg;
		std::string _strIndex;

	public:
		static const uint8_t charshift;
		static const uint8_t gapChar;

		/**
		 * Dimer constructor: states c1 and c2 separated by d wildcards.
		 */
		UniversalKmer(const int d, const uint8_t c1, const uint8_t c2);

		/**
		 * Trimer constructor
		 */
		UniversalKmer(const int g1, const int g2, const uint8_t c1, const uint8_t c2, const uint8_t c3);

		/**
		 * Construct from any AbstractKmer
		 */
		UniversalKmer(const AbstractKmer&);

		void mutate(const int offset, const uint8_t state);

		uint8_t charAt(const int n) const{ return cg[n].chr; }

		int gapsAfter(const int n) const{ return cg[n].gap;	}

		motif_columns_type getMotifColumns() const{
			motif_columns_type motif_columns;
			int pos = 0;
			motif_columns.push_back(pos++);
			for(int i=1; i<numMatches(); i++){
				pos += gapsAfterDirect(i-1);
				motif_columns.push_back(pos++);
			}
			return motif_columns;
		}

		/* functions duplicates that get inlined when base class is known */
		uint8_t charAtDirect(const int n) const{return cg[n].chr;}
		int gapsAfterDirect(const int n) const{return cg[n].gap;}

		UniversalKmer* clone() const {
			return new UniversalKmer(*this);
		}

		std::string getStrIndex() const {
			return _strIndex;
		}

		std::string showStrIndex() const {
			std::stringstream s;
			for (size_t i=0; i<_strIndex.length(); ++i) {
				s << "[" << (int)_strIndex[i] << "]";
			}
			return s.str();
		}

};

std::ostream& operator<<(std::ostream& os, const UniversalKmer &k);

#endif /* UNIVERSALKMER_H_ */
