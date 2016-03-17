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

#include "AlignStream.h"

#ifndef __CLUSTALW_H__
#define __CLUSTALW_H__

class ClustalW: public AlignStream {
private:

public:
	ClustalW() :
		AlignStream() {
	}
	ClustalW(const AlignStream& a) :
		AlignStream(a) {
	}
	~ClustalW() {
	}
	bool read(std::istream & s, bool checkFormat);
	std::ostream & print(std::ostream &s) const;
	ClustalW& operator=(const AlignStream& rhs) {
		this->assign(rhs.begin(), rhs.end());
		return *this;
	}
};

#endif
