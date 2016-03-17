#ifndef STATELIB_H_
#define STATELIB_H_

#include <functional>
#include <limits>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;

#include "alphabet.h"
#include "columnState.h"
#include "utils.h"

class StateLib {
	public:
		typedef std::list<std::pair<double, uint8_t> > pair_list_type;
		typedef std::vector<pair_list_type> sigStates_type;
	private:
		typedef std::vector<ColumnState> state_vec_type;
		state_vec_type states;
		Alphabet alphabet;

		/**
		 * threshold for a state considered to be relevant for a certain aa,
		 * i.e., if score(aa, state) >= thresh
		 */
		double significanceThreshold;

		/**
		 * type of column state
		 */
		StateType type;

		static const Alphabet dummyAlphabet;
	public:

		StateLib() : alphabet(dummyAlphabet), significanceThreshold(std::numeric_limits<double>::infinity()), type(UNDEF) {}
		StateLib(const Alphabet &A) :
			alphabet(A), significanceThreshold(std::numeric_limits<double>::infinity()), type(UNDEF) {
		}

		static void initStates(const StateLib &states, const StateType &type);

		size_t size() const {
			return states.size();
		}

		const ColumnState& operator[](const int i) const {
			return states[i];
		}

		ColumnState& operator[](const int i) {
			return states[i];
		}

		void push_back(const ColumnState &s) {
			states.push_back(s);
		}

		state_vec_type::const_iterator begin() const {
			return states.begin();
		}

		state_vec_type::const_iterator end() const {
			return states.end();
		}

		void pop_back() {
			states.pop_back();
		}

		double getSignificanceThreshold() const {
			return significanceThreshold;
		}

		StateType getStateType() const {
			return type;
		}

		void setType(const StateType &t) {
			type = t;
		}

		/**
		 * A vector containing, for each amino acid, a list of pairs
		 * (score, state index), sorted in descending order.
		 * This can be used to quickly find all states with at least
		 * score x for amino acid a.
		 */
		sigStates_type sigStates;
		void computeSigStates(const double thresh);

		void sortByEntropy();

		void writeTikz(const int per_line, const std::string filename) const;

};

#endif /* STATELIB_H_ */
