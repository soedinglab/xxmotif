/*
 * alphabet.h
 *
 *  Created on: Oct 12, 2009
 *      Author: eckhart
 */

#ifndef MPRO_ALPHABET_H_
#define MPRO_ALPHABET_H_

#include "../alphabet.h"

#include <string>
#include <vector>

#include <iostream>
using std::cerr;
using std::endl;

/**
 * Representation of an alphabet of characters.
 */
class Alphabet
{

public:

	typedef char charType;

	static const charType startStopChar;

    Alphabet()
    {}

    /**
     * Construct alphabet from Holger's C-structure.
     */
    Alphabet(const a_type &A);
    const int& size() const
    {
        return numCharacters;
    }

    const charType& wildcardChar() const
    {
        return alphabet[wildcardIndex()];
    }

    const int& wildcardIndex() const
    {
        return _wildcardIndex;
    }

    const charType& itoc(const int i) const
    {
        return alphabet[i];
    }

    const int& ctoi(const charType c) const
    {
        return let2code[c];
    }

    const charType& operator[](const int i) const
    {
        return alphabet[i];
    }

private:
    int _wildcardIndex;
    std::vector<char> alphabet;
    std::vector<int> let2code;
    int numCharacters;
};

#endif /* ALPHABET_H_ */
