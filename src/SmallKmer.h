#ifndef SMALLKMER_H_
#define SMALLKMER_H_


#include "alphabet.h"

#include "AbstractKmer.h"

#include <stdint.h>
#include <string>
#include <vector>

class SmallKmer: public AbstractKmer {
	private:
		static std::string convBase(uint64_t v, long base);
		static int chrBits;
		static int gapBits;
		static const uint64_t allMask;
		static std::vector<uint64_t> chrMask;
		static std::vector<uint64_t> gapMask;
		static std::vector<uint64_t> prefixMask;
		static std::vector<uint64_t> suffixMask;
		static std::vector<int> bitsBeforeChar;
		static std::vector<int> bitsBeforeGap;
		static bool initialized;

	public:
		static int maxNumMatches;
		static void init(int cbits, int gbits);

		SmallKmer(const int gap, const uint8_t c1, const uint8_t c2);

		SmallKmer(const int gap1, const int gap2, const uint8_t c1, const uint8_t c2, const uint8_t c3);

		void mutate(const int offset, const uint8_t state);

		uint8_t charAt(const int n) const{
			assert(n<numMatches());
			return (uint8_t)((id.val.num & chrMask[n]) >> bitsBeforeChar[n]);
		}

		int gapsAfter(const int n) const{
			assert(n>=0);
			return n+1 < numMatches() ? static_cast<int>((id.val.num & gapMask[n]) >> bitsBeforeGap[n]) : 0;
		}

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
		uint8_t charAtDirect(const int n) const{
			return static_cast<uint8_t>((id.val.num & chrMask[n]) >> bitsBeforeChar[n]);
		}
		int gapsAfterDirect(const int n) const{
			return static_cast<int>((id.val.num & gapMask[n]) >> bitsBeforeGap[n]);
		}

		SmallKmer* clone() const{
			return new SmallKmer(*this);
		}

		std::string asBinaryString() const;
};

#endif /* SMALLKMER_H_ */
