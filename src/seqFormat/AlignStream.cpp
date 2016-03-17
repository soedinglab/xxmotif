/* Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

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

#ifndef __ALIGNSTREAM_INL_H__
#define __ALIGNSTREAM_INL_H__

#include "AlignStream.h"

AlignStream::AlignStream(const std::vector<FormatContainer> & _data) {
	data.assign(_data.begin(), _data.end());
}

AlignStream::~AlignStream()
{}

void AlignStream::assign(const_iterator beg, const_iterator end) throw (SeqException)
{
	if(beg == end) throw(SeqException("no valid sequences found"));
	data.assign(beg, end);

	size_t len = beg->GetSeq().length();

	for (AlignStream::const_iterator itr = beg + 1; itr != end; ++itr) {
		if (itr->GetSeq().length() != len) {
			throw(SeqException("data elements have different lengths or different identifiers"));
		}
	}
}

std::istream & AlignStream::ReadThroughLine(std::istream & s) {
	char ch;
	while (s.get(ch)) {
		if (ch == '\n')
			return s;
	}
	return s;
}


#endif
