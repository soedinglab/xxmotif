#ifndef _EXTENSION_H
#define _EXTENSION_H

#include <list>

#include <unordered_map>

#include "Kmer.h"
#include "../NullModel.h"
#include "../alphabet.h"


typedef std::list<std::shared_ptr<Kmer>, Pool_alloc<std::shared_ptr<Kmer> > > elongList;
typedef std::unordered_map<AbstractKmer::id_type, float, std::hash<AbstractKmer::id_type> > fullHashType;

class Extension {

public:

	/**
	 * Options to control the visualization output of the motif extension phase
	 */
	typedef struct {
		std::ostream *stream;
		bool showEdgesToLosers; /* whether to show edges to nodes other than the best ones */
		bool showEdgesToLookups; /* whether to visualize lookups in hash table of already extended motifs */
		std::string id_suffix; /* string representation of the root motif, needed for unique ids in .dot file */
	} graphVizOptions;

	typedef std::unordered_map<AbstractKmer::id_type, std::shared_ptr<Kmer>, std::hash<AbstractKmer::id_type> > hashmap;

	static void initialize(const MProGlobal::SeqInfo_t &S);
	static void initialize();

	static void calibrateMotif(Kmer &kres, bool debug = false);

	static MProGlobal::SeqInfo_t S;

	static elongList getBestExtensions(
			const elongList &krs,
			fullHashType *extendedMotifs = NULL,
			const graphVizOptions* vizopts = NULL, bool debugOut = false);

	static const char nucleotide[5];
	static const uint8_t consBases[11][5];

private:
	template <typename KmerType>
	static double getKmerProb(KmerType& Kmer);

	static void sumKmerProbRec(const unsigned char* kmer, unsigned char* expanded_kmer, int pos, const int length, double& p);
};

template <typename KmerType>
double Extension::getKmerProb(KmerType& Kmer){
	static unsigned char expanded_kmer[100];
	static unsigned char kmer[100];

	int i=0, j=0;
	for(; i<Kmer.numMatches()-1; i++){
		kmer[j++] = Kmer.charAt(i);
		for(int gaps = 0; gaps < Kmer.gapsAfter(i); gaps++){ kmer[j++] = 10; }
	}
	kmer[j] = Kmer.charAt(i);

	/* loop through all possible nucleotides for first IUPAC character */
	int pos = 0;
	unsigned char charAt = kmer[pos];
	int nbDegenerations = consBases[charAt][0];

	double p = 0;
	for(int loop = 1; loop <= nbDegenerations; loop++){
		expanded_kmer[pos] = consBases[charAt][loop];
		sumKmerProbRec(kmer, expanded_kmer, pos + 1, Kmer.length(), p);
	}
	return p;
}

#endif
