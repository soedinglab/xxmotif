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

#include <map>
#include "CustomMultFasta.h"
#include <iterator>
#include <algorithm>
#include <sstream>

void CustMultFasta::read(std::istream & s)
{
	char ch;

	if (!(s >> ch)) {
		s.setstate(std::ios::badbit);
		return;
	} else{
		s.putback(ch);
	}

	s >> ch;

	if (ch != '>') {
		throw badFormat("file not in customized mulitple FASTA format: Header does not start with \">\"");
	}

	std::string id;
	std::string seq;
	while (1) {
		s.get(ch);
		if (ch == '\n')
			break;
		id += ch;
	}

	seq.reserve(1000);
	std::vector<FormatContainer> _data;
	std::string temp;

	int nb = 0;
	while (!s.eof()) {
		s >> ch;
		if (ch == '>') {
			s.putback(ch);
			break;
		}else {
			seq += ch;
			std::getline(s, temp);
			seq += temp;
		}

		if(seq.at(seq.length()-1) == '$'){
			seq.erase(seq.length()-1);
			_data.push_back(FormatContainer(id, seq));
			std::ostringstream stm;
			stm << ++nb;
			id = stm.str();
			seq.clear();
		}
	}

	this->assign(_data.begin(), _data.end());
}

std::ostream & CustMultFasta::print(std::ostream &s) const {
	CustMultFasta::const_iterator i = this->begin(), j = this->end();
	size_t len = i->GetSeq().length(), k = 0;
	s << "CLUSTAL W" << "\n\n";
	while (k < len) {
		size_t offset = (k + 60 < len) ? k + 60 : k + (len - k);
		for (i = this->begin(); i < j; ++i) {
			s << i->GetId() << '\t';
			std::copy(i->GetSeq().begin() + k, i->GetSeq().begin() + offset,
					std::ostream_iterator<char>(s, ""));
			s << '\n';
		}
		s << '\n';
		k = offset;
	}
	return s;
}
