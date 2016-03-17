/*
 * utils.h
 *
 *  Created on: Oct 14, 2009
 *      Author: eckhart
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "alphabet.h"
#include "sequence.h"
#include "../utils.h"

#include <algorithm>
#include <list>
#include <string>
#include <valarray>
#include <vector>
#include <stdint.h>

/**
 * For a filename given, return its basename, i.e. the prefix up to
 * the last '.' character.
 */
std::string getBasename(const char* f);

template<class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
	os << "[";
	for (size_t i = 0; i < v.size(); ++i) {
		os << v[i];
		if (i < v.size() - 1) {
			os << ", ";
		} else {
			os << "]";
		}
	}
	return os;
}


template<class T>
std::ostream& operator<<(std::ostream &os, const std::valarray<T> &v) {
	os << "[";
	for (size_t i = 0; i < v.size(); ++i) {
		os << v[i];
		if (i < v.size() - 1) {
			os << ", ";
		} else {
			os << "]";
		}
	}
	return os;
}

template<class T, class Alloc>
std::ostream& operator<<(std::ostream &os, const std::list<T, Alloc> &l) {
	os << "[";
	size_t count = 1;
	for (typename std::list<T, Alloc>::const_iterator it = l.begin(); it != l.end(); ++it, ++count) {
		os << *it;
		if (count < l.size()) {
			os << ", ";
		}
	}
	os << "]";
	return os;
}

template<class T, class S>
std::ostream& operator<<(std::ostream &os, const std::pair<T, S> &p) {
	os << "(" << p.first << "," << p.second << ")";
	return os;
}

template<typename T, typename S, typename CmpT = std::less<T>,
		typename CmpS = std::less<S> >
class pairCmp {
public:
	bool operator()(const std::pair<T, S> &p1, const std::pair<T, S> &p2) {
		if (CmpT()(p1.first, p2.first)) {
			return true;
		} else {
			if (CmpT()(p2.first, p1.first)) {
				return false;
			} else {
				return (CmpS()(p1.second, p2.second));
			}
		}
	}
};

template<typename I> bool isSortedContainer(I begin, I end) {
	I prev = begin;
	I curr = prev;
	++curr;
	while (curr != end) {
		if (*curr < *prev) {
			return false;
		}
		++curr;
		++prev;
	}
	return true;
}

template<class T, class C> int removeDuplicatesFromSortedContainer(T &cont, C eq) {
	int num_dups = 0;
	typename T::iterator current = cont.begin();
	while (current != cont.end()) {
		typename T::iterator next = current;
		++next;
		while (next != cont.end() && eq(*current, *next)) {
			next = cont.erase(next);
			++num_dups;
		}
		++current;
	}
	return num_dups;
}

#endif /* UTILS_H_ */
