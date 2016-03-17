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

#ifndef __ALIGNSTREAM_H__
#define __ALIGNSTREAM_H__

#include "SeqExceptions.h"
#include <utility>
#include <string>
#include <vector>
#include "FormatContainer.h"

class AlignStream {
private:
	typedef std::pair<std::string, std::string> baseType;
	std::vector<FormatContainer> data;

protected:
	std::istream & ReadThroughLine(std::istream &);

public:
	AlignStream(const std::vector<FormatContainer> & _data);
	AlignStream(const AlignStream &a) { this->assign(a.begin(), a.end()); }
	AlignStream(void) {}
	virtual ~AlignStream();

	int size(void) const { return static_cast<int>(data.size()); }
	FormatContainer& operator[](const int i) { return data[i];	}
	const FormatContainer& operator[](const int i) const { 	return data[i];	}

	typedef std::vector<FormatContainer>::iterator iterator;
	typedef std::vector<FormatContainer>::const_iterator const_iterator;
	iterator begin() { return data.begin();	}
	iterator end()   { return data.end();   }
	const_iterator begin() const { return data.begin();	}
	const_iterator end()   const { return data.end();	}

	const std::vector<FormatContainer>& Data(void) { return data;}

	void assign(const_iterator beg, const_iterator end) throw (SeqException);
};

#endif
