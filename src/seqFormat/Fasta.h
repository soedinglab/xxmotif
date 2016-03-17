/*

 Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

 Remove the brackets to email me.

 This file is part of libsequence.

 libsequence is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 libsequence is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

 */

/*!
 \class Sequence::Fasta Sequence/Fasta.hpp
 \ingroup seqio
 Publicly derived from Sequence::Seq, this class defines
 how to read and print sequences in FASTA format, which looks like:\n
 >sequence name 1\n
 ATGATGATCAGATAGACATAGCAGATACATGT\n
 >sequence name 2\n
 ATGTTGGTTTTTTTTTAGAGATGTTTATAGGT\n
 ETC...
 
 @short FASTA sequence stream
 */

#ifndef FASTA_H
#define FASTA_H

#include "FormatContainer.h"
#include "../alphabet.h"

class Fasta : public FormatContainer {
private:
	const a_type A;

public:
	Fasta(const a_type alph) : A(alph) {}
	Fasta(const std::string &id, const std::string &seq, const a_type alph);
	Fasta(const char *id, const char *seq, const a_type alph);
	~Fasta(){}

	bool read(std::istream &s) throw (badFormat, std::exception);
	std::ostream& print(std::ostream& s) const;
};
#endif
