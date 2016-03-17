#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <iosfwd>
#include <string>
#include <utility>
#include "Alignment.h"
#include "../alphabet.h"

class FormatContainer{
public:
	FormatContainer(){}
	FormatContainer(const std::string& id, const std::string& seq) : _id(id), _seq(seq) {}
	FormatContainer(const char* id, const char* seq) : _id(id), _seq(seq){}
	FormatContainer(const FormatContainer& fc) : _id(fc.GetId()), _seq(fc.GetSeq()){}
	virtual ~FormatContainer() {}

	std::string GetId() const { return _id; }
	std::string GetSeq() const { return _seq; }

	int length() const {return static_cast<int>(_seq.length());}

protected:
	std::string _id;
	std::string _seq;

private:

};

#endif
