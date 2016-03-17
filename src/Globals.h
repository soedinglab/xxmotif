#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "alphabet.h"
#include "seqset.h"
#include "log.h"
#include "ThresholdChecker.h"

#include "getopt_pp/getopt_pp.h"
#include "memoryPool/pool_alloc.h"
#include <string.h>
#include <stdlib.h>
#include <unordered_set>
#include <math.h>
#include <sys/stat.h>


#ifdef NDEBUG
const bool kDebug = false;
#else
const bool kDebug = true;
#endif  // DEBUG

#define MAX_CONS_MUTATIONS 30			/* maximum number of mutations calibrated in conservation track */
#define MAX_MUTATED_FRACTION 1  		/* maximum fraction of mutated residues */
#define MAX_CONSERVATION_LENGTH 15		/* maximum calibrated length of conservation track				*/
#define PWM_LENGTH 30

static const int mutIndex[4] = {(int)pow(2,0), (int)pow(2,4), (int)pow(2,8), (int)pow(2,12)}; /* indices for A, C, G, T mutations in conservation track */

enum motif_type{
	ALL,
	PALINDROME,
	TANDEM,
	NOPALINDROME,
	NOTANDEM,
	FIVEMERS,
	NO_VALID_MOTIF_TYPE
};

enum merge_type{
	LOW,
	MEDIUM,
	HIGH,
	NO_VALID_MERGE_MODE
};
std::ostream& operator<<(std::ostream &os, const merge_type &m);

enum TerminusMode {
	NONE, POS, NEG, BOTH
};

/* which kind of supplementary information to use */
enum SuppInfMode {
	SUPP_NO, SUPP_DISOCONS, SUPP_NNET
};
std::ostream& operator<<(std::ostream &os, const SuppInfMode &v);

enum Types_t {
	CS_63, BLOSUM45_21, BLOSUM62_21, BLOSUM80_21, GONNET_21, UNDEF
};

class StateType {
public:
	Types_t type;
	StateType() : type(UNDEF) {}
	StateType(const Types_t &t) : type(t) {}
	bool operator==(const Types_t &t) const {
		return type == t;
	}
	bool operator!=(const Types_t &t) const {
		return !(*this==t);
	}
	static std::string toString(const StateType &t) {
		switch (t.type) {
		case CS_63 : return std::string("CS_63");
		case BLOSUM45_21 : return std::string("BLOSUM45_21");
		case BLOSUM62_21 : return std::string("BLOSUM62_21");
		case BLOSUM80_21 : return std::string("BLOSUM80_21");
		case GONNET_21 : return std::string("GONNET_21");
		case UNDEF : return std::string("UNDEF");
		default : std::cerr << "Unknown StateType!" << std::endl; exit(1);
		}
	}
	static Types_t fromString(const std::string &s) {
		if (s == "CS_63") return CS_63;
		else if (s=="BLOSUM45_21") return BLOSUM45_21;
		else if (s=="BLOSUM62_21") return BLOSUM62_21;
		else if (s=="BLOSUM80_21") return BLOSUM80_21;
		else if (s=="GONNET_21") return GONNET_21;
		else { assert(false); return UNDEF; }
	}
	static int sizeOfAlphabet(const StateType &t) {
		switch (t.type) {
		case CS_63 : return 63;
		case BLOSUM45_21 : return 21;
		case BLOSUM62_21 : return 21;
		case BLOSUM80_21 : return 21;
		case GONNET_21 : return 21;
		case UNDEF : return std::numeric_limits<int>::min();
		default : std::cerr << "Unknown StateType!" << std::endl; exit(1);
		}
	}
};

typedef std::list<int, Pool_alloc<int> > motif_columns_type;

class Global{
public:
	Global(int nopt, char *options[]);
	~Global();

	static a_type 		A;				/** Alphabet used **/
	static ss_type 		posSet;			/** positive sequence set **/
	static ss_type 		negSet;			/** negative sequence set **/

	static bool 		usePositionalProbs;
	static int 			backgroundOrder;

	static bool			positionalProbsRanking;

	static int			GAPS;			/** number of gap combinations allowed for the start Motifs **/
	static merge_type	mergeMode;
	static double		gapOpening;
	static double		gapExtension;
	static int 			maxMultipleSequences;
	static int 			conservationLength;
	static int			maxPosSetSize;
	static double*		motifNbCorrection;
	static ThresholdChecker instanceThreshold;	/* threshold for final decision whether an instance is a motif or not */
	static bool         removeHomology; /* remove perfectly conserved motifs with maxMatchPostions */

	static double		consCorrection; 	  /* pValCorrection for conservation pValue */
	static double		overrepCorrection; 	  /* pValCorrection for overrepresentatino pValue */
	static double		consPvalWeight;		  /* weight for pVal combination */
	static int			maxSeqCount;

	static int 			maxMotifsPerSequence;	  /* maximum number of motifs per sequence */

	static bool 	    useRankPvalues;
	static bool			useAliFree;

	static bool			maximizeMotifLength;
	static bool			noRefinementPhase;

	static bool   		multipleOccurrence; /** multiple instances should be found in one sequence **/
	static bool   		oneOccurrence; /** multiple instances should be found in one sequence **/
	static bool			zeroOrOneOccurrence;
	static bool			ungappedOutput; /** output file in pwm-folder has gapped output **/
	static bool			revcomp;		/** search on reverse complement of sequences **/
	static bool 		repeatFiltering;
	static bool			lowComplexityFilter;

	/****** startMotif ****/
	static char* 		startMotif;		/** start motif for motif discovery **/
	static char*		profFile;		/** file with a the start profile for motif discovery **/
	static int			startRegion;	/** start of enriched region **/
	static int 			endRegion;		/** end of enriched region **/
	static motif_type 	type;			/** type of motif: normal, palin or sequential **/
	static seq_format	seqFormat;		/** format of input sequence: FASTA, CLUSTALW **/

  	static float***		conservationProbs;	/** conservation probability for a nmer with mutations in n species **/
  	static float***		alignmentFreeProbs;	/** conservation probability for a nmer with mutations in n species **/

  	static double* 		posBg_log;		/* logarithm of distribution of positive set **/
  	static double* 		posBg;			/* background of positive set **/
  	static double* 		negBg_log;		/* logarithm of distribution of negative set **/
  	static double*		negBg;			/* background of negative set **/

    static double		pseudo;			/** pseudocounts for pwm **/
    static double		plusFrac;		/** plus fraction in motif iteration **/

    static int			neff_pwm;			/** effective number of different bases in one pwm column **/
    static int 			neff_discrete; /** effective alphabet size before switching to PWM */

	static int 			downstream;		/** distance between the alignment point and the end of the input sequences **/

	static char* 		outputDirectory; /** output Directory **/
	static char*		name;			/** input file name **/
	static char*		shortFileName;	/** input file name without path and .**/
	static char*		negFile;		/** negative input file name **/
	static char* 		benchmarkFolder;/** folder in which the benchmark results are stored **/
	static char* 		pwmFolder;		/** folder in which pwms for bulyk benchmark are stored **/

	static int			maxMotifLevel; /** max number of extensions accepted for next level */
	static double		minCoverage; /** min number of sequences with instance of a motif */

	static int			minMatchPositions; /** max number of non-wildcard motif positions */
	static int			maxMatchPositions; /** min number of non-wildcard motif positions */

	/* background model order >= 0 */
	static int order;
	/* background model pseudocounts factor >= 0 */
	static float pseudocountsFactor;
	/* background model counts offset >= 0 */
	static float countsOffset;

	/*
	 * Expectiation Maximization (EM) parameters
	 * */

	/*
	 * higher-order modeling and EM mode
	 */
	static bool	em;

	/*
	 * Parameters to initialize models from file
	 * */

	/*
	 * binding sites to initialize a single Markov model. Sequence lengths must
	 * not differ and be provided line-by-line
	 */
	static char* bindingSiteFile;
	/*
	 * binding sites lengths
	 */
	static int bindingSiteLength;

	/*
	 * Markov model file name (without ending) to initialize a single Markov
	 * model. Files .conds and .probs need to be available
	 */
	static char* markovModelFile;
	/*
	 * Markov model size
	 */
	static int markovModelLength;

	/*
	 * Parameters to initialize models from XXmotif results. Parameters
	 * bindingSiteFile and markovModelFile must not be provided simultaneously
	 * */

	/*
	 * number of one or more XXmotif models in the ranking used to initialize
	 * Markov models. The remaining parameters avaiable to choose models from
	 * XXmotif results are ignored
	 */
	static std::vector<int> nrModels;
	/*
	 * min. number of XXmotif models used to initialize Markov models.
	 * Independent on pValueThreshold and minOccurrence
	 */
	static int minModels;
	/*
	 * max. number of XXmotif models used to initialize Markov models
	 */
	static int maxModels;
	/*
	 * max. p-value of XXmotif models used to initialize Markov models. Not
	 * applied to min. number of models
	 */
	static double pValueThreshold;
	/*
	 * min. percentage of sequences containing a binding site instance. Not
	 * applied to min. number of models
	 */
	static float minOccurrence;
	/*
	 * add columns to the left and right of XXmotif models used to initialize
	 * Markov models
	 */
	static std::vector<int> addColumns;
	/*
	 * use model-specific specificity factor by calculating the percentage of
	 * positive sequences containing a corresponding binding site instance
	 */
	static bool msq;

	/*
	 * Markov model parameters
	 * */

	/*
	 * Markov model order
	 */
	static int modelOrder;
	/*
	 * Markov model pseudo-counts factor(s). Markov model order k fixes vector size
	 * to k+1
	 */
	static std::vector<float> alpha;

	/*
	 * learn hyper parameter alpha during interpolation of markov models
	 */
    static bool learnHyperParameter;

    /*
     * flag in order to print output for debugging he alpha learning part only
     */
    static bool debugAlphalearning;

    /*
     *  whether or not to use position specific alhpa's
     */
    static bool positionSpecificAlphas;


	/*
	 * interpolate between higher- and lower-order probabilities
	 */
	static bool interpolate;

	/*
	 * Interpolated Markov background model parameters
	 * */

	/*
	 * Background model order
	 */
	static int modelOrderBg;
	/*
	 * Background model pseudo-counts factor
	 */
	static float alphaBg;

	/*
	 * EM parameters
	 * */

	/*
	 * initialize Markov model but skip EM phase
	 */
	static bool noExpectationMaximizationPhase;
	/*
	 * specificity factor approximates the percentage of sequences contributing to
	 * the Markov model
	 */
	static float q;
	static float qmax;
	/*
	 * EM convergence parameter
	 */
	static float epsilon;
	/*
	 * likelihood or max. order model parameter EM convergence
	 */
	static bool likelihoodConvergence;
	/*
	 * max. number of EM iterations
	 */
	static int maxEMIterations;
	/*
	 * update interpolated Markov model probabilities with last EM iteration's
	 * pseudo-counts
	 */
	static bool lastCondsPseudoCounts;
	/*
	 * calculate 0th-order interpolated Markov model pseudo-counts from initial
	 * 0th-order probabilities
	 */
	static bool monoProbsPseudoCounts;
	/*
	 * calculate 0th-order interpolated Markov model pseudo-counts from initial
	 * 0th-order probabilities using 0th-order pseudo-counts factor N * q *
	 * alphaZeroFactor
	 */
	static float alphaZeroFactor;

	/*
	 * Weighting parameters
	 */

	/*
	 * intensity or significance values for positive sequences. The higher the
	 * values the higher the weights
	 */
	static char* sequenceIntsFile;
	/*
	 * parameter to initialize models from XXmotif results by weighting
	 * instances with corresponding sequence weigths. Option --sequenceIntsFile
	 * must be provided simultaneously. Options --bindingSiteFile and
	 * --markovModelFile must not be provided simultaneously
	 */
	static bool initInts;
	/*
	 * rank-based weighting. Defaults to intensity-based weighting
	 */
	static bool rankWeighting;
	/*
	 * quantile to estimate the background intensity value (or rank). Sequences
	 * having their intensity value (rank) below (above) the background
	 * intensity value (rank) get assigned to weight zero. Defaults to 0
	 */
	static float backgroundQuantile;
	/*
	 * background intensity value. Sequences having their intensity value below
	 * the background intensity value get assigned to weight zero. Defaults to
	 * the min. intensity value
	 */
	static float backgroundIntensity;
	/*
	 * background intensity rank. Sequences having their intensity rank above
	 * the background intensity rank get assigned to weight zero. Defaults to
	 * the max. rank
	 */
	static float backgroundRank;

	/*
	 * Binding site weighting parameters
	 */

	/*
	 * intensity or significance values for binding site sequences. The higher
	 * the values the higher the weights. Parameter bindingSiteFile must be
	 * provided simultaneously
	 */
	static char* bindingSiteIntsFile;
	/*
	 * binding site rank-based weighting. Defaults to intensity-based weighting
	 */
	static bool bindingSiteRankWeighting;
	/*
	 * quantile to estimate the background intensity value (or rank). Binding
	 * sites having their intensity value (rank) below (above) the background
	 * intensity value (rank) get assigned to weight zero. Defaults to 0
	 */
	static float bindingSiteBackgroundQuantile;
	/*
	 * background intensity value. Binding sites having their intensity value
	 * below the background intensity value get assigned to weight zero.
	 * Defaults to the min. binding site intensity value
	 */
	static float bindingSiteBackgroundIntensity;
	/*
	 * background intensity rank. Binding sites having their intensity rank
	 * above the background intensity rank get assigned to weight zero. Defaults
	 * to the max. binding site rank
	 */
	static float bindingSiteBackgroundRank;

	/*
	 * Scoring parameters
	 * */

	/*
	 * evaluate model(s) on training sequences
	 */
	static bool testPosSequences;
	/*
	 * evaluate model(s) on background sequences
	 */
	static bool testNegSequences;
	/*
	 * evaluate model(s) on sequences in FASTA format. Specify one or more files.
	 * Sequence lengths may differ
	 */
	static std::vector<std::string> testSequenceFile;
	/*
	 * evaluate PWM model(s) used to initialize Markov model(s) on test
	 * sequences
	 */
	static bool evaluatePWMs;

	/*
	 * calculate log probabilities instead of log likelihood ratios
	 */
	static bool logProbs;

	/*
	 * Output parameters
	 * */

	/*
	 * save Markov models after initialization to file
	 */
	static bool saveInitModels;
	/*
	 * save Markov models after EM phase to file
	 */
	static bool saveModels;
	/*
	 * save EM iteration's sequence likelihoods and positional odds to file
	 */
	static bool saveExpectationMaximizationLikelihoods;
	/*
	 * save EM iteration's Markov models to file
	 */
	static bool saveExpectationMaximizationModels;
	/*
	 * verbose printouts
	 */
	static bool	verbose;

	/*
	 * Internal parameters
	 */

	/*
	 * monomer background frequencies
	 */
	static double* freqs;
	/*
	 * calculate background probabilities for k-mers with gaps. Gaps are mandatory
	 * in order to initialize from XXmotif's models
	 */
	static bool gaps;
	/*
	 * list of sequences to evaluate model(s)
	 */
	static std::list<ss_type> testSet;

	/* csblast stuff */
	static int			cswlen;
	static std::string	csprofiles;
	static int 			csbest;

	/* MadonaPro stuff */
	static double		disoconsWeight; /** weight for P value combination */
	static TerminusMode termMode; /** append dollar to pos_set sequence termini? */
	static float		dollarScaleFactor; /** factor to scale $ frequencies */
	static StateType	type_of_states;		/** type of profile states */
	static int			maxIterations;	/** max number of pwm iterations */
	typedef std::unordered_set<std::string> tracked_t;
	static tracked_t trackedMotifs;	/** all dimer/trimer substrings will be tracked, i.e. [KR][DE]EL */
	static bool isTracked(const std::string &s);
	static bool 		trackedElongation;
	static bool			trackedOnly;	/** use tracked motifs only as initial seeds */
	static double		extensionECut; /** initial seeds with E<=cut will be extended */
	static int			extensionMinCut; /** minimum number of extended seeds (if enough) */
	static int			extensionMaxCut; /** maximum number of extended seeds */
	static double		aaStateSigThresh; /** min score s[a] to consider a state relevant for amino acid s */
	static double		aaSeqFreqThresh; /** min fraction of sequences with aa conserved for extension */
	static bool			batch_mode; /** if running non-interactively (suppress progress indicators) */

	static bool fixedPosition;
	static double finalFilterThreshold; /** E-value threshold for filtering final motifs */

	static bool 		DEBUG;
	static std::string	argv0; /** name of executable as given on the command line */

	static std::string nnetFilename;

	static char* String(const char *s);
	static void createDirectory(const char *s);
private:
	bool readCommandLineOptions(int argc, char *argv[]);
	void printHelpOutput();
};

inline char* Global::String(const char *s){
  return strdup(s);
}

inline void Global::createDirectory(const char *s){
	struct stat St;
	if( stat( s, &St ) != 0 ){
		fprintf(stderr, "output directory does not exist\n");
		char* command = (char*)calloc(1024, sizeof(char));
		sprintf(command, "mkdir %s", s);
		if(system(command) != 0){
			fprintf(stderr, "Directory %s could not be created\n", s);
			exit(-1);
		}
		free(command);
	}
}

namespace GetOpt
{
	template <> inline _Option::Result convert<char*>(const std::string& s, char*& d, std::ios::fmtflags)
	{
		_Option::Result ret = _Option::BadType;
		d = Global::String(s.c_str());
		ret = _Option::OK;
		return ret;
	}

	template <> inline _Option::Result convert<motif_type>(const std::string& s, motif_type& d, std::ios::fmtflags)
	{
		_Option::Result ret = _Option::BadType;

		if(s.compare("PALINDROME") == 0) {
			d = PALINDROME; std::cerr << "Motif Type: PALINDROME" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("TANDEM") == 0) {
			d = TANDEM; std::cerr << "Motif Type: TANDEM" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("NOPALINDROME") == 0){
			d = NOPALINDROME; std::cerr << "Motif Type: NOPALINDROME" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("NOTANDEM") == 0){
			d = NOTANDEM; std::cerr << "Motif Type: NOTANDEM" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("ALL") == 0){
			d = ALL; std::cerr << "Motif Type: ALL" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("FIVEMERS") == 0){
			d = FIVEMERS; std::cerr << "Motif Type: FIVEMERS" << std::endl;
		}else {
			d = NO_VALID_MOTIF_TYPE; std::cerr << "\nERROR: Motif type \"" << s << "\" not possible" << std::endl;
		}
		return ret;
	}

	template <> inline _Option::Result convert<merge_type>(const std::string& s, merge_type& d, std::ios::fmtflags)
	{
		_Option::Result ret = _Option::BadType;

		if(s.compare("LOW") == 0) {
			d = LOW; std::cerr << "Similarity threshold for merging motifs set to: LOW" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("MEDIUM") == 0) {
			d = MEDIUM; std::cerr << "Similarity threshold for merging motifs set to: MEDIUM" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("HIGH") == 0){
			d = HIGH; std::cerr << "Similarity threshold for merging motifs set to: HIGH" << std::endl;
			ret = _Option::OK;
		}else {
			d = NO_VALID_MERGE_MODE; std::cerr << "\nERROR: Similarity threshold for merging motifs \"" << s << "\" not possible" << std::endl;
		}
		return ret;
	}
	template <> inline _Option::Result convert<seq_format>(const std::string& s, seq_format& d, std::ios::fmtflags)
	{
		_Option::Result ret = _Option::BadType;

		if(s.compare("CLUSTALW") == 0) {
			d = CLUSTALW; std::cout << "Format of Input Sequences: CLUSTALW" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("FASTA") == 0){
			d = FASTA; std::cout << "Format of Input Sequences: FASTA" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("MFASTA") == 0){
			d = MFASTA; std::cout << "Format of Input Sequences: multiple FASTA" << std::endl;
			ret = _Option::OK;
		}else if(s.compare("CUSTOM") == 0){
			d = CUSTOM; std::cout << "Format of Input Sequences: customized multiple FASTA" << std::endl;
			ret = _Option::OK;
		}else {
			d = NO_VALID_SEQ_FORMAT; std::cerr << "\nERROR: Sequence Format \"" << s << "\" not accepted" << std::endl;
		}
		return ret;
	}
}

inline bool Global::isTracked(const std::string &s) {
   for (auto m = trackedMotifs.begin(); m != trackedMotifs.end(); ++m) {
      if (s.compare(*m) == 0) {
        return true;
      }
  }
  return false;
}

#endif /* _GLOBALS_H_ */
