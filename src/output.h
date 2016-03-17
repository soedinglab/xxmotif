#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "refinementPhase/MotifContainer.h"

class Output{
public:
	static void printTopN(const MotifContainer& motifs, int N = 10);
	static void printMergedMotifs(const MotifContainer& motifs, int N = 10);
	static void printFinalMotifs(const MotifContainer& motifs, int N = 10);

	static void printOutput(const MotifContainer& motifs, int N = 10);
	static void printBenchmark(const MotifContainer& motifs);
	static void print_pwm_folder(const MotifContainer& motifs);
	static void print_cs_output(const MotifContainer& motifs);
	static void writeMTFFile(const MotifContainer& motifs);
	static void writeRFile(const MotifContainer& motifs);

	static void writeBestKmerSites( MotifContainer& motifs, int minSites, char* baseFileName );
private:
	static void writeBlocksFile(const MotifContainer& motifs);
	static void printMotifDistribution(const MotifContainer& motifs);
	static void writeSummaryFile(const MotifContainer& motifs);
	static void printStartPosFile(const MotifContainer& motifs);
	static void printPvalFile(const MotifContainer& motifs);
	static void writeMemeFile(const MotifContainer& motifs);
	static void printSequenceFile();
	static void printPosSetFreq();	
	static void printScoreDistribution(const MotifContainer& motifs);
	static void printDistanceDistribution(const MotifContainer& motifs);
};


#endif /* OUTPUT_H_ */
