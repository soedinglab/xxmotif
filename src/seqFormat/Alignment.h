#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <utility>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "SeqExceptions.h"
#include "../alphabet.h"

enum seq_format{
	FASTA,
	CLUSTALW,
	MFASTA,
	CUSTOM,
	NO_VALID_SEQ_FORMAT
};

typedef struct {
	int n; 		/* length of entity sequence */
	unsigned short mseq; 	/* number of sequences in multiple alignment */
	unsigned char **S; 		/* multiple sequence alignment*/
	char **info; 			/* description of entity */
	float weight; 			/* weight of entity */
} sequence_type;
typedef sequence_type *e_type;

typedef struct {
	char		*name;		/* input filename for entities */
	e_type		*entity;	/* array of sequence entities */
	int			nent;		/* number of input entities */
	int	        max_MultSeq;/* maximum number of sequences in multiple alignment */
	int			max_leng;	/* sequence entity maximumlength */
	int			min_leng;	/* sequence entity minimum length */
	double*		avgLength;  /* avgLength of sequence set */
	int*		total;		/* total # residues */
} seqset_type;
typedef seqset_type	*ss_type;

struct SequenceData{
	SequenceData(std::string ID, std::string Seq) : id(ID), seq(Seq){}
	std::string id;
	std::string seq;
};

typedef std::list<SequenceData> alignmentType;

class Alignment {

public:
	Alignment(char* filename, a_type A, seq_format format, const bool addTermini, int max_seqs);

	ss_type getSeqSet(){ return P; }

private:
	void readFasta(std::ifstream &inFile);
	void readClustalW(std::ifstream &inFile);
	void readMultipleFasta(std::ifstream &inFile);
	void readCustomizedFasta(std::ifstream &inFile);

	void pushBackAlig(alignmentType& alig);
	
	e_type getEntity(const alignmentType& align) const;

	char* filename;
	const a_type A;
	
	ss_type P;
	int maxSeqAlloc;
	const bool addTermini;
	int max_seqs;
};

#endif
