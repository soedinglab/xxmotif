#include "Fasta.h"
#include <stdexcept>
#include <iostream>
#include <string>

Fasta::Fasta(const std::string &id, const std::string &seq, const a_type alph) :
	FormatContainer(id, seq), A(alph){
}

Fasta::Fasta(const char *id, const char *seq, const a_type alph) : FormatContainer(id, seq), A(alph){
}

bool Fasta::read(std::istream & stream) throw (badFormat,	std::exception) {
	bool newBlock = false;
	char ch;
	bool seqflag;

	if (!(stream >> ch)) {
		stream.setstate(std::ios::badbit);
		return static_cast<bool>(stream);
	} else
		stream.putback(ch);

	std::string temp;

	stream >> ch;

	if (ch != '>') {
		throw badFormat("file not in FASTA format: Header does not start with \">\" : " + ch);
	}
	_id.clear();
	while (1) {
		stream.get(ch);
		if (ch == '\n')
			break;
		_id += ch;
	}

	seqflag = 1;
	_seq.clear();
	_seq.reserve(1000);
	while (seqflag) {
		stream >> ch;
		if (ch == '/'){
			std::getline(stream, temp);
			seqflag = 0;
			newBlock = true;
		}else if (ch == '>') {
			stream.putback(ch);
			seqflag = 0;
		} else if (stream.eof())
			seqflag = 0;
		else {
			_seq += ch;
			std::getline(stream, temp);

			if(temp.length() > 0 && temp.at(temp.length()-1) == '$'){
				std::string alphstr(A->alphabet, A->n + 1);
				std::cerr << alphstr << std::endl;
				if (alphstr.find('$') == std::string::npos) {
					throw badFormat("file not in FASTA format: sequence ends with $ => use --format CUSTOM");
				}
			}
			_seq += temp;
		}
	}
	return newBlock;
}

std::ostream & Fasta::print(std::ostream & stream) const {
	stream << '>' << _id << '\n' << _seq;
	return stream;
}
