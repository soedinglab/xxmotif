#include <iostream>
#include <string>
#include "./AbstractKmer.h"
#include "Globals.h"
#include "SmallKmer.h"
#include "elongationPhase/Kmer.h"
#include "elongationPhase/Match.h"
#include "elongationPhase/elongationCandidates.h"
#include "memoryPool/pool_alloc.h"

using std::cerr;
using std::endl;

std::string AbstractKmer::toString(const int numBestCols, const char* alphabet) const {
  std::string s;
  if (alphabet == NULL) {
    for (uint8_t i = 0; i < numMatches(); ++i) {
      s += MProGlobal::getS().states[charAt(i)].toString(numBestCols);
      if (i + 1 < numMatches()) {
        s += std::string(gapsAfter(i),
            MProGlobal::getS().alphabet.wildcardChar());
      }
    }
  } else {
    for (uint8_t i = 0; i < numMatches(); ++i) {
      s += alphabet[charAt(i)];
      if (i + 1 < numMatches()) {
        s += std::string(gapsAfter(i), '.');
      }
    }
  }
  return s;
}

std::string AbstractKmer::bestAAString() const {
  std::string s(length(), MProGlobal::getS().alphabet.wildcardChar());
  int strPos = 0;
  int matchPos = 0;
  while (matchPos < numMatches()) {
    s[strPos] = MProGlobal::getS().states[charAt(matchPos)].nthSignificantChar(
        0);
    if (matchPos + 1 < numMatches()) {
      strPos += 1 + gapsAfter(matchPos);
    }
    ++matchPos;
  }
  return s;
}

std::string AbstractKmer::bestNucString() const {
  std::string s(length(), 'N');
  int strPos = 0;
  int matchPos = 0;
  while (matchPos < numMatches()) {
    s[strPos] = ElongCandidates::IUPAC_CHARS[charAt(matchPos)];
    if (matchPos + 1 < numMatches()) {
      strPos += 1 + gapsAfter(matchPos);
    }
    ++matchPos;
  }
  return s;
}

bool AbstractKmer::isMatchPosition(const int n) const {
  if (n < 0 || n >= length()) {
    return false;
  } else {
    int pos = 0;
    for (uint8_t i = 0; i + 1 < numMatches(); ++i) {
      if (pos == n) {
        return true;
      } else if (pos > n) {
        return false;
      }
      pos += 1 + gapsAfter(i);
    }
    return pos == n;
  }
}

bool AbstractKmer::operator==(const AbstractKmer &other) const {
  if (numMatches() != other.numMatches() || length() != other.length()) {
    return false;
  } else {
    for (int i = 0; i < numMatches(); ++i) {
      if (charAt(i) != other.charAt(i)) {
        return false;
      }
    }
    for (int i = 0; i + 1 < numMatches(); ++i) {
      if (gapsAfter(i) != other.gapsAfter(i)) {
        return false;
      }
    }
  }
  return true;
}

std::ostream& operator<<(std::ostream& os, const AbstractKmer &k) {
  os << k.toString();
  return os;
}

std::ostream& operator<<(std::ostream& os, const AbstractKmer::id_type &id) {
  if (id.isNumeric) {
    os << id.val.num;
  } else {
    os << std::string(id.val.str);
  }
  return os;
}
