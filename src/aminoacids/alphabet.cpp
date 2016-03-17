#include "alphabet.h"

const char Alphabet::startStopChar = '$';

Alphabet::Alphabet(const a_type &A) :
	_wildcardIndex(0), alphabet(&A->alphabet[0], &A->alphabet[A->n + 1]),
			let2code(), numCharacters(A->n + 1) {

	/*
	 * Everything can be copied upon initialization,
	 * except for let2code, here we have to determine the highest
	 * possible index first.
	 */
	char maxChar = 0;
	for (int i = 0; i < numCharacters; ++i) {
		if (alphabet[i] > maxChar) {
			maxChar = alphabet[i];
		}
	}
	let2code.resize(maxChar + 1);
	for (int i = 0; i < numCharacters; ++i) {
		let2code[alphabet[i]] = i;
	}

}
