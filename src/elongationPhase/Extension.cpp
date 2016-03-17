#include <unordered_set>

#include "Extension.h"
#include "ExtensionTable.h"
#include "elongationCandidates.h"
#include "../AbstractKmer.h"
#include "../UniversalKmer.h"
#include "../SmallKmer.h"
#include "../pValCalculation.h"
#include "../prodpval_stat.h"
#include "../AbstractKmer.h"
#include "../Globals.h"
#include "../SmallKmer.h"
#include "../elongationPhase/Kmer.h"
#include "../elongationPhase/Match.h"
#include "../elongationPhase/elongationCandidates.h"
#include "../memoryPool/pool_alloc.h"

using std::vector;
using std::string;

const char Extension::nucleotide[5] = {'A','C','G','T','N'};

const uint8_t Extension::consBases[11][5] = {
	{1, 1,0,0,0}, /* A: loop = 1, return A */
	{1, 2,0,0,0}, /* C: loop = 1, return C */
	{1, 3,0,0,0}, /* G: loop = 1, return G */
	{1, 4,0,0,0}, /* T: loop = 1, return T */
	{2, 1,2,0,0}, /* M: loop = 2, return A/C */
	{2, 1,3,0,0}, /* R: loop = 2, return A/G */
	{2, 1,4,0,0}, /* W: loop = 2, return A/T */
	{2, 2,3,0,0}, /* S: loop = 2, return C/G */
	{2, 2,4,0,0}, /* Y: loop = 2, return C/T */
	{2, 3,4,0,0}, /* K: loop = 2, return G/T */
	{1, 0,0,0,0}  /* N: loop = 4, return N */
};

MProGlobal::SeqInfo_t Extension::S;

void Extension::initialize() {
		S.D_MAX = 4;
		S.alphabet = Alphabet(Global::A);
		S.states = StateLib(S.alphabet);

		/* create states */
		S.states.push_back(ColumnState(1, 0, 0, 0)); /* A */
		S.states.push_back(ColumnState(0, 1, 0, 0)); /* C */
		S.states.push_back(ColumnState(0, 0, 1, 0)); /* G */
		S.states.push_back(ColumnState(0, 0, 0, 1)); /* T */
		S.states.push_back(ColumnState(0.5, 0.5, 0, 0)); /* M */
		S.states.push_back(ColumnState(0.5, 0, 0.5, 0)); /* R */
		S.states.push_back(ColumnState(0.5, 0, 0, 0.5)); /* W */
		S.states.push_back(ColumnState(0, 0.5, 0.5, 0)); /* S */
		S.states.push_back(ColumnState(0, 0.5, 0, 0.5)); /* Y */
		S.states.push_back(ColumnState(0, 0, 0.5, 0.5)); /* K */

		vector<uint8_t> v;
		S._posSet.push_back(v); /* first sequence is at index 1 */
		for (int i = 1; i <= ElongCandidates::nbSequences; i++) {
			e_type &ent = Global::posSet->entity[i];
			vector<uint8_t> v((uint8_t*) (ent->S[0]), (uint8_t*) (ent->S[0]
					+ ent->n + 1));
			S._posSet.push_back(Sequence(v));
		}
		ExtensionTable::initialize(S, log(0.4), 0.2, 0.4); // all states are used (0.5 > 0.4)
}

elongList Extension::getBestExtensions(
		const elongList &krs,
		fullHashType *extendedMotifs,
		const graphVizOptions *vizopts, bool debugOut) {

	int firstCharSeq = 1;

	static int hashHitCount;
	static int calibratedCount;

	/* remember which motifs for the next stage have already been calibrated */
	std::unordered_set<AbstractKmer::id_type, std::hash<AbstractKmer::id_type> > calibratedMotifs;
	//calibratedMotifs.rehash(1e4);

	/* list of candidates for the next stage */
	elongList nextCandidates;

	ExtensionTable::stateList_type estates;

//	cerr << "\nBEG********************************\nelongList\n";
//	for (elongList::const_iterator kres = krs.begin(); kres != krs.end(); ++kres) {
//		cerr << **kres << endl;
//	}
//	cerr << "END********************************\n\n";

	for (elongList::const_iterator kres = krs.begin(); kres != krs.end(); ++kres) {
		/* defines whether current motif is only in list because it was found in hash over all calibrated motifs */
		bool hashMotif = false;
		if((*kres)->seeds.begin() == (*kres)->seeds.end()){
			hashMotif = true;
		}
		int candidates = 0;
		if (vizopts && vizopts->stream) {
			*(vizopts->stream) << '\t' << '"'
					<< (*kres)->getKmer()->getId().toString()
					<< vizopts->id_suffix << '"' << " [shape=box,label=" << '"';
				*(vizopts->stream) << (*kres)->getKmer()->toString(0,
						ElongCandidates::IUPAC_CHARS) << "\\n";
			*(vizopts->stream) << exp((*kres)->p_set);
			if(!hashMotif){ *(vizopts->stream) << "\\n" << (*kres)->setSize;}
			*(vizopts->stream)  << '"' << "];" << std::endl << std::flush;

		}
		if (debugOut) {
				cout << "Extending " << (*kres)->getKmer()->toString(0,
						ElongCandidates::IUPAC_CHARS) << endl << std::flush;
		}

		const int kmerLength = (*kres)->getKmer()->length();

		/*create list with possible extension positions sorted by there significance */
		std::list<int> interestingPositions;

		const int testPositions = kmerLength-1;

		for(int i=1; i<testPositions; i++){
			if ((*kres)->getKmer()->isMatchPosition(i)) continue;
			interestingPositions.push_back(i);
		}
		for(int i=1; i<= S.D_MAX; i++){
			for(int j=-1; j<=1; j+=2) interestingPositions.push_back(i*j > 0 ? i*j + kmerLength -1 : i*j);
		}

		for(std::list<int>::iterator it = interestingPositions.begin(); it != interestingPositions.end(); it++){
			int32_t offset = static_cast<int32_t>(*it);

			if (debugOut){
				cout << " at offset " << offset << ": ";
				cout	<< (*kres)->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS) << endl << std::flush;
			}

			estates.clear();
			if(!hashMotif){
				/* reset extension table and count amino acids in the matches at the new position*/
				ExtensionTable::reset();

				MatchContainer::const_iterator mit = (*kres)->seeds.begin();
				int num_seeds = 0;
				while (mit != (*kres)->seeds.end()) {
					if(mit->pos < (*kres)->enrichment.startRegion || mit->pos > (*kres)->enrichment.endRegion){
						// match not in enriched region
					}else if (mit->pos + offset >= (int) S._posSet[mit->seq].size() || mit->pos + offset < 0) {
						// match protrudes sequence boundaries, will be eliminated later
					}else {
						ExtensionTable::count(mit->seq, S._posSet[mit->seq][mit->pos + offset], Global::multipleOccurrence);
						++num_seeds;
					}
					++mit;
				}
				/* would never be significant if only one seed is left */
				if (num_seeds < 2) {
					continue;
				}
				/* get promising extensions and construct new Kmer for calibration from each one */
				int setSize = (*kres)->setSize;
				ExtensionTable::getExtensionStates(setSize, estates, debugOut);
			}else{
				/* find possible extensions by checking all possible extension in hash */
				for (uint8_t s = 0; s < S.states.size(); ++s) {
					AbstractKmer* tmpKmer;
					if ((*kres)->getKmer()->numMatches() == SmallKmer::maxNumMatches) {
						tmpKmer = new UniversalKmer(*(*kres)->getKmer());
					}else{
						tmpKmer = (*kres)->getKmer()->clone();
					}
					tmpKmer->mutate(offset, s);
//					cerr << "searching in hash for " << tmpKmer->toString() << " ... " << tmpKmer->getId() << endl;
//					cerr << (long)extendedMotifs << endl;
//					if (extendedMotifs==NULL) {
//						cerr << "ERROR: extendedMotifs is NullPointer" << endl;
//						exit(1);
//					}
//					cerr << "hash contains " << extendedMotifs->size() << " elements: " << endl;
//					for (fullHashType::const_iterator dit = extendedMotifs->begin(); dit != extendedMotifs->end(); ++dit) {
//						cerr << *dit << ", ";
//					}
					fullHashType::iterator it = extendedMotifs->find(tmpKmer->getId());
					if( it != extendedMotifs->end()){
						estates.push_back(s);
					}
					delete tmpKmer;
				}
			}

			/* if no states can be elongated, or motif has too many significant states for branch and bound */
			if (!estates.empty() && (*kres)->getKmer()->numMatches() < Global::maxMatchPositions) {
				if (debugOut) {
					size_t count = 0;
					cout << "possible extensions: [";
					for (ExtensionTable::stateList_type::const_iterator
							it = estates.begin(); it != estates.end(); ++it) {
							cout << ElongCandidates::IUPAC_CHARS[(int) *it];
						++count;
						if (count < estates.size()) cout << ", ";
					}
					cout << "]" << endl;
				}

				/* try extensions, calibrate motif and add to list of candidates if better than ancestor */
				for (ExtensionTable::stateList_type::const_iterator
						state_it = estates.begin(); state_it != estates.end(); ++state_it) {
					std::shared_ptr<Kmer> cand(new Kmer(*(*kres)->getKmer()));
					AbstractKmer* kmer = cand->getKmer();
					kmer->mutate(offset, *state_it);

					/* if longer than branch and bound can take, continue */
					if(kmer->length() > Global::maxMatchPositions){continue;}

					if (debugOut){
						cout << "\t" << cand->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS) << ": " << std::flush;
					}

					bool hashHit = false;
					if(extendedMotifs){
						AbstractKmer::id_type id = kmer->getId();
						fullHashType::iterator it = extendedMotifs->find(id);
						/* motif is already calibrated */
						if( it != extendedMotifs->end()){
							hashHit = true;
							/* motif has not already been seen in this round */
							if(calibratedMotifs.find(id) == calibratedMotifs.end()){
								double pVal = it->second;
								if (pVal <= (*kres)->p_set){
									cand->p_set = pVal;
									nextCandidates.push_back(cand);
								}else if (vizopts && vizopts->stream && hashMotif) {
									*(vizopts->stream) << '\t' << '"'
									<< cand->getKmer()->getId().toString()
									<< vizopts->id_suffix << '"'
									<< " [shape=box,label=" << '"';
										*(vizopts->stream) << cand->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS);
									*(vizopts->stream) << "\\n" << exp(pVal);
									*(vizopts->stream)  << '"' << ", color=\"#999999\", fontcolor=\"#999999\"];" << std::endl << std::flush;
								}
								calibratedMotifs.insert(id);
								//cerr << "!!! found already in another round: " << cand->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS) << ": " << exp(pVal) << endl;
							}
						}
					}else if(calibratedMotifs.find(kmer->getId()) != calibratedMotifs.end()) {
						/* motif has already been calibrated, do nothing */
						hashHitCount++;
						hashHit = true;
					}
					cand->enrichment = (*kres)->enrichment;

					int32_t moveStartPos = (offset > 0 ? 0 : offset);
					cand->enrichment.startRegion += moveStartPos;
					cand->enrichment.endRegion += moveStartPos;
					cand->enrichment.max += moveStartPos;

					if(hashHit){
						if (vizopts && vizopts->stream
						&& vizopts->showEdgesToLookups) {
							*(vizopts->stream) << '\t' << '"'
									<< (*kres)->getKmer()->getId().toString()
									<< vizopts->id_suffix << '"' << " -> "
									<< '"'
									<< cand->getKmer()->getId().toString()
									<< vizopts->id_suffix << '"'
									<< " [style=dotted]" << std::endl << std::flush;
						}
						if (debugOut) cout << "found in hash " << endl	<< std::flush;
					}else{
						MatchContainer::iterator mit = (*kres)->seeds.begin();
						while (mit != (*kres)->seeds.end()) {
							int32_t newStart = static_cast<int32_t>(mit->pos + moveStartPos);
							int32_t newEnd = static_cast<int32_t>(mit->pos + offset);

							if (newStart >= firstCharSeq && newEnd < (int)S._posSet[mit->seq].size()) {
								// mit->pos must be the first character of the mutated matched kmer
								const int base = S._posSet[mit->seq][mit->pos + offset];
									switch(*state_it){
									/*M*/ case 4: if(base == 1 || base == 2) cand->seeds.push_back(Match(mit->seq, newStart));break;
									/*R*/ case 5: if(base == 1 || base == 3) cand->seeds.push_back(Match(mit->seq, newStart));break;
									/*W*/ case 6: if(base == 1 || base == 4) cand->seeds.push_back(Match(mit->seq, newStart));break;
									/*S*/ case 7: if(base == 2 || base == 3) cand->seeds.push_back(Match(mit->seq, newStart));break;
									/*Y*/ case 8: if(base == 2 || base == 4) cand->seeds.push_back(Match(mit->seq, newStart));break;
									/*K*/ case 9: if(base == 3 || base == 4) cand->seeds.push_back(Match(mit->seq, newStart));break;
							  /*A,C,G,T*/ default:if(base == *state_it+1) cand->seeds.push_back(Match(mit->seq, newStart));break;
									}
							}
							++mit;
						}

						//if (debugOut) cerr << "New candidate (before cali): " << *cand << endl;
						calibrateMotif(*cand, debugOut);
						//if (debugOut) cerr << "New candidate (calibrated):  " << *cand << endl;

						calibratedCount++;
						if (cand->p_set <= (*kres)->p_set){
							candidates++;
							if (debugOut)
									cout << "ACCEPTED: " << cand->p_set
									<< " <= " << (*kres)->p_set << endl
									<< std::flush;
							nextCandidates.push_back(cand);
							if (vizopts && vizopts->stream
									&& vizopts->showEdgesToLookups) {
								*(vizopts->stream) << '\t' << '"'
										<< (*kres)->getKmer()->getId().toString()
										<< vizopts->id_suffix << '"' << " -> "
										<< '"'
										<< cand->getKmer()->getId().toString()
										<< vizopts->id_suffix << '"'
										<< std::endl << std::flush;
							}
						} else {
							if (debugOut) cout << "rejected" << endl
									<< std::flush;
							if (vizopts && vizopts->stream
									&& vizopts->showEdgesToLosers) {
								*(vizopts->stream) << '\t' << '"'
										<< (*kres)->getKmer()->getId().toString()
										<< vizopts->id_suffix << '"' << " -> "
										<< '"'
										<< cand->getKmer()->getId().toString()
										<< vizopts->id_suffix << '"'
										<< " [style=dashed, color=\"#999999\"]"
										<< std::endl;
								*(vizopts->stream) << '\t' << '"'
										<< cand->getKmer()->getId().toString()
										<< vizopts->id_suffix << '"'
										<< " [shape=box,label=" << '"';
									*(vizopts->stream) << cand->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS);
								*(vizopts->stream) << "\\n" << exp(cand->p_set);
								if(!hashMotif) *(vizopts->stream) << "\\n" << cand->setSize;
								*(vizopts->stream)  << '"' << ", color=\"#999999\", fontcolor=\"#999999\"];" << std::endl << std::flush;
							}
						}

						if(extendedMotifs)(*extendedMotifs)[kmer->getId()] = static_cast<float>(cand->p_set);
						calibratedMotifs.insert(kmer->getId());
					}
				}
			} else{
				if (debugOut) cout << endl << std::flush;
			}
			if(offset >= kmerLength && candidates > 0){
				break;
			}
		}
	}

	elongList bestUnextended;
	if (!krs.empty()) {
		elongList::const_iterator k_it = krs.begin();
		double best_p = (*k_it)->p_set;
		while (k_it != krs.end() && (*k_it)->p_set == best_p) {
			bestUnextended.push_back(*k_it);
			++k_it;
		}
	}

	if (!nextCandidates.empty()) {
		/* sort candidates, remove duplicates and retain only the best N ones */
		nextCandidates.sort(Kmer::lessPtr);
		elongList::iterator cand_it = nextCandidates.begin();

		int count = 0;
		while (count < Global::maxMotifLevel && cand_it != nextCandidates.end()) {
			elongList::iterator next = cand_it;
			++next;
			while (next != nextCandidates.end() && (**cand_it).isDuplicate(**next)) {
				next = nextCandidates.erase(next);
			}
			++count;
			++cand_it;
		}
		while (cand_it != nextCandidates.end()) {
			if (vizopts && vizopts->stream) {
				*(vizopts->stream) << '\t' << '"'
						<< (*cand_it)->getKmer()->getId().toString()
						<< vizopts->id_suffix << '"'
						<< " [shape=plaintext,style=filled,label=" << '"';
					*(vizopts->stream) << (*cand_it)->getKmer()->toString(0, ElongCandidates::IUPAC_CHARS);
				*(vizopts->stream) << "\\n" << exp((*cand_it)->p_set);
				if((*cand_it)->seeds.begin() != (*cand_it)->seeds.end()) *(vizopts->stream) << "\\n" << (*cand_it)->setSize;
				*(vizopts->stream)  << '"' << ", color=\"#EEEEEE\", fontcolor=\"#000000\"];" << std::endl << std::flush;

			}
			cand_it = nextCandidates.erase(cand_it);
		}
		if (debugOut){
			cout << "Calibrated: " << calibratedCount << "   Hash hits: "
				 << hashHitCount << endl << std::flush;
			cout << "Total: " << nextCandidates.size()
				 << " candidates for next stage" << endl << endl << std::flush;
		}
		elongList extensions =	getBestExtensions(nextCandidates, extendedMotifs, vizopts,debugOut);
		if (!extensions.empty()) {
			if (bestUnextended.empty() || extensions.front()->p_set	<= bestUnextended.front()->p_set) {
				return extensions;
			}
		}
	}
	return bestUnextended;
}

void Extension::calibrateMotif(Kmer &kres, bool debug) {
	//cout << "number Of Seeds: " << kres.seeds.size() << endl;

	AbstractKmer& kmer = *(kres.getKmer());
	double p = 1;
	if(kmer.numMatches() > SmallKmer::maxNumMatches){
		UniversalKmer* uk = reinterpret_cast<UniversalKmer*>(&kmer);
		p = getKmerProb<UniversalKmer>(*uk);
	}else{
		SmallKmer* sk = reinterpret_cast<SmallKmer*>(&kmer);
		p = getKmerProb<SmallKmer>(*sk);
	}

	//fprintf(stderr, "%s: %d\n", kres.getKmer()->toString(0, ElongCandidates::IUPAC_CHARS).c_str(), kres.getKmer()->numMatches());
	p *= pow(Global::overrepCorrection, kres.getKmer()->numMatches());
	kres.p_pos = p;

	/* pVal for finding the motif in a random sequence (use average length of sequences with motif) */
	int startRange = kres.enrichment.startRegion;
	int endRange = kres.enrichment.endRegion;

	int kmerLength = kres.getKmer()->length();

	int maxSeq = ElongCandidates::nbSequences;
	double LOG_Bonferonni = calculate_log_bonferonni(kmer.getMotifColumns(), LogTable::LOG_i[Global::neff_discrete]);

	int count = 0;
	int N = ElongCandidates::nbSequences;
	memset(ElongCandidates::motifsPerSequenceCount, 0, Global::posSet->nent*sizeof(int));
	if(Global::multipleOccurrence){
		if(!Global::revcomp){
			/* do not count overlapping motifs */
			if(Global::usePositionalProbs){
				N = (endRange - startRange + 1) * ElongCandidates::nbSequences;

				MatchContainer& startPositions = kres.seeds;
				int oldSeq = 0; int oldPos = 0;
				for(MatchContainer::iterator it = startPositions.begin(); it != startPositions.end(); it++){
					if(oldSeq == it->seq && oldPos + kmerLength > it->pos) continue;
					if(it->pos >= startRange && it->pos <= endRange){
						ElongCandidates::motifsPerSequenceCount[it->seq]++;
						if(ElongCandidates::motifsPerSequenceCount[it->seq]<= Global::maxMotifsPerSequence){
							count++;
						}
						oldSeq = it->seq; oldPos = it->pos;
					}
				}
			}else{
				N = ElongCandidates::totalPositions - ElongCandidates::nbSequences*kmerLength;
				MatchContainer& startPositions = kres.seeds;
				int oldSeq = 0; int oldPos = 0;
				for(MatchContainer::iterator it = startPositions.begin(); it != startPositions.end(); it++){
					if(oldSeq == it->seq && oldPos + kmerLength > it->pos) continue;
					ElongCandidates::motifsPerSequenceCount[it->seq]++;
					if(ElongCandidates::motifsPerSequenceCount[it->seq]<= Global::maxMotifsPerSequence){
						count++;
					}
					oldSeq = it->seq; oldPos = it->pos;
				}
			}
		}else{
			N = ElongCandidates::totalPositions - ElongCandidates::nbSequences*kmerLength;
			MatchContainer& startPositions = kres.seeds;
			int oldSeq = 0; int oldPos = 0;
			ss_type set = Global::posSet;
			//if(debug)fprintf(stderr, "\n\n");
			for(MatchContainer::iterator it = startPositions.begin(); it != startPositions.end(); it++){
				/* don't count overlapping motifs on the same strand */
				if(oldSeq == it->seq && oldPos + kmerLength > it->pos) continue;

				int seqLength = set->entity[it->seq]->n;
				if(Global::revcomp && it->pos > seqLength / 2){
					int revPos = seqLength - it->pos - kmerLength + 2;
					//if(debug)fprintf(stderr, "\n (%d/%d), revPos: %d\t", it->seq, it->pos, revPos);
					MatchContainer::iterator it_rev = it;
					bool found = false;
					while(!found && it_rev != startPositions.begin()){
						it_rev--;
						if(it_rev->seq != it->seq || it_rev->pos <= revPos + kmerLength - 1){
							found = true;
						}
					}
					//if(debug)fprintf(stderr, "newRev: %d/%d, compare: %d/%d\n", it_rev->seq, it_rev->pos, it->seq, revPos-kmerLength+1);
					/* don't count overlapping motifs on the reverse strands */
					if(found && it_rev->seq == it->seq && it_rev->pos > revPos - kmerLength + 1){
					   // if(debug) fprintf(stderr, " => skip\n");
						 continue;
					}
				}
				ElongCandidates::motifsPerSequenceCount[it->seq]++;
				if(ElongCandidates::motifsPerSequenceCount[it->seq]<= Global::maxMotifsPerSequence){
					count++;
				}
				//if(debug)fprintf(stderr, "%d/%d\t", it->seq, it->pos);
				oldSeq = it->seq; oldPos = it->pos;
			}
		}
		//if(debug)fprintf(stderr, "\n\n");
	}else{
		if(Global::usePositionalProbs) p = 1-pow(1-p, endRange-startRange+1);
		else p = 1-pow(1-p, ElongCandidates::avgLength);

		bool* seqOcc = (bool*)calloc(maxSeq+1, sizeof(bool));
		MatchContainer& startPositions = kres.seeds;
		for(MatchContainer::iterator it = startPositions.begin(); it != startPositions.end(); it++){
			if(!seqOcc[it->seq] && it->pos >= startRange && it->pos <= endRange){
				seqOcc[it->seq] = true;
				count++;
			}
		}
		free(seqOcc);
	}
	kres.setSize = count;
	double pVal = calculateOrderStatisticsPvalue(count, N, p);

	if(Global::posSet->max_MultSeq > 1){
		PVal_Calculator& pValCalculator = PVal_Calculator::getInstance();
		const motif_columns_type& columns = kmer.getMotifColumns();
		double log_prod_pcons = 0;

		const MatchContainer& seeds = kres.seeds;
		double pConsCorrection = pValCalculator.getConsCorrection(static_cast<int>(columns.size()));
		for(MatchContainer::const_iterator it = seeds.begin(); it != seeds.end(); it++){
			double pValCons = pValCalculator.calculatePvalCons(it->seq, it->pos, 0, columns);
			pValCons *= pConsCorrection;
			if(pValCons > 1){
				pValCons = 1;
			}else{
				//if(debug)fprintf(stderr, "based: %f, free: %f\n", pValAliBased, pValAliFree);
				log_prod_pcons += log(pValCons);
			}
		}
		double pValOverRep = pVal;
		double pValCons = pValCalculator.getCombinedPval_log(log_prod_pcons, static_cast<int>(seeds.size()));
		/* combined pValue */
		pVal = pValCalculator.getWeightedPval_log(pValOverRep, pValCons, Global::consPvalWeight);

		if(debug){
			fprintf(stderr, "\n%s: pVal: %e, pValCons: %e => %e\n", kmer.bestNucString().c_str(), exp(pValOverRep), exp(pValCons), exp(pVal));
		}
	}

	kres.p_set = pVal + LOG_Bonferonni;

	if(debug)fprintf(stderr, "\ncount: %d, N: %d, p_pos: %e, pVal: %e, Bonferonni: %e, p_set: %e\n", kres.setSize, N, kres.p_pos, exp(pVal), exp(LOG_Bonferonni), exp(kres.p_set));
}

void Extension::sumKmerProbRec(const unsigned char* kmer, unsigned char* expanded_kmer, int pos, const int length, double& p){
	if(pos == length){
		p += NullModel::getProbability(expanded_kmer, length);
	}else{
		unsigned char charAt = kmer[pos];
		int nbDegenerations = consBases[charAt][0];

		for(int loop = 1; loop <= nbDegenerations; loop++){
			expanded_kmer[pos] = consBases[charAt][loop];
			sumKmerProbRec(kmer, expanded_kmer, pos + 1, length, p);
		}
	}
}
