#include "Alignment.h"
#include "Fasta.h"
#include "Clustalw.h"
#include "CustomMultFasta.h"
#include <string.h>
#include "../Globals.h"

using namespace std;

Alignment::Alignment(char* file, a_type alphabet, seq_format format, const bool addTermini, int max_seqs) :
		filename(file), A(alphabet), addTermini(addTermini), max_seqs(max_seqs) {
	
	ifstream inFile;
	inFile.open(filename);
	
	P = (ss_type)malloc(1*sizeof(seqset_type));
	P->name = (char *)malloc((strlen(filename)+2)*sizeof(char)); 
	strcpy(P->name,filename);
	
	maxSeqAlloc = 20000; // is reallocated if not sufficient
	P->entity = (e_type*) malloc((maxSeqAlloc+1)*sizeof(e_type));
	P->nent=0;
	P->min_leng = numeric_limits<int>::max();
	P->max_leng = 0;
	P->max_MultSeq = 0;
	
	if (!inFile) {
	    cout << "Unable to open file " << filename << endl;
	    exit(1);
	}
	if(format == FASTA){
		readFasta(inFile);
	}else if(format == CLUSTALW){
		readClustalW(inFile);
	}else if(format == MFASTA){
		readMultipleFasta(inFile);
	}else if(format == CUSTOM){
		readCustomizedFasta(inFile);
	}
	inFile.close();
	
	P->total = (int*)malloc((P->nent+1)*sizeof(int));
	P->avgLength = (double*)malloc((P->nent+1)*sizeof(double));
}

void Alignment::pushBackAlig(alignmentType& alig){
	P->nent++;
	if(P->nent  > maxSeqAlloc){
		maxSeqAlloc += 20000;
		P->entity = (e_type*) realloc(P->entity, (maxSeqAlloc+1)*sizeof(e_type));
	}
	e_type E = getEntity(alig);
	P->entity[P->nent] = E;
	if(E->n > P->max_leng){ P->max_leng = E->n; }
	if(E->n < P->min_leng){ P->min_leng = E->n; }
	if(E->mseq > P->max_MultSeq){ P->max_MultSeq = E->mseq; }
	
	alig.clear();
}

void Alignment::readFasta(ifstream &inFile){
	while (!inFile.fail()) {
		if(P->nent == max_seqs){ break; }
		Fasta x(A);
		try{
			if(x.read(inFile))
				throw badFormat("input file not in FASTA format: \n\tif you have a file with several alignments seperated by // use option --format MFASTA");
		}catch( SeqException str ){
	        cout << "\nError:\t" << str << "\n\n"; exit(-1);
	    }

        alignmentType alig;
        alig.push_back(SequenceData(x.GetId(), x.GetSeq()));

        pushBackAlig(alig);
	}
}

void Alignment::readMultipleFasta(ifstream &inFile){
	int nent = 0;
	bool newBlock = true;
	alignmentType alig;
	while (!inFile.fail()) {
		Fasta x(A);
		try{
			newBlock = x.read(inFile);
		}catch( SeqException str ){
	        cout << "Exception raised in Alignment Block " << nent << ": " << str << '\n'; exit(-1);
	    }

		if(x.GetSeq().size() == 0)	continue;

		alig.push_back(SequenceData(x.GetId(), x.GetSeq()));

		if(newBlock == true || inFile.fail()){
			pushBackAlig(alig);
			if(++nent == max_seqs){ break; }
		}
	}
}

void Alignment::readClustalW(ifstream &inFile){
	int nent = 0;
	bool checkFormat = true;
	while (!inFile.fail()) {
		ClustalW x;
		try{
			checkFormat = x.read(inFile, checkFormat);
		}catch( SeqException str ){
			cout << "Exception raised in Alignment Block " << nent << ": " << str << '\n';
			exit(-1);
		}

		alignmentType alig;
		for (int i = 0 ; i < x.size() ; i++){
			if(x[i].GetSeq().size() == 0) continue;
			alig.push_back(SequenceData(x[i].GetId(), x[i].GetSeq()));
		}

		pushBackAlig(alig);
		if(++nent == max_seqs){ break; }
	}
}

void Alignment::readCustomizedFasta(ifstream &inFile){
	int nent = 0;

	while (!inFile.fail()) {
		CustMultFasta x;
		try{
			x.read(inFile);
		}catch( SeqException str ){
			cout << "Exception raised in Alignment Block " << nent << ": " << str << '\n';
			exit(-1);
		}

		alignmentType alig;
		for (int i = 0 ; i < x.size() ; i++){
			if(x[i].GetSeq().size() == 0) continue;
			alig.push_back(SequenceData(x[i].GetId(), x[i].GetSeq()));
		}
		pushBackAlig(alig);
		if(++nent == max_seqs){ break; }
	}
}



e_type Alignment::getEntity(const alignmentType& align) const {
	/* allocate memory for entity */
	e_type E = (e_type) malloc(1 * sizeof(sequence_type));

	E->n = static_cast<int>(align.front().seq.length());
	if(addTermini){ E->n += 2; }

	E->mseq = static_cast<unsigned short>(align.size());

	E->info = (char **) malloc(E->mseq * sizeof(char*));
	E->S = (unsigned char**) malloc(E->mseq * sizeof(char*));

	int m=0;
	for(alignmentType::const_iterator it = align.begin(); it != align.end() && m < E->mseq; it++, m++){
		E->info[m] = (char*)calloc((it->id.length()+1), sizeof(char));
		E->S[m] = (unsigned char*)calloc(E->n+1, sizeof(unsigned char));

		string id = it->id;
		/* remove carriage return at end of sequence if existing */
		const size_t pos = id.find_last_of("\r");
		if(pos != std::string::npos) id.erase(pos);

		/* copy sequence id */
		strcpy(E->info[m], id.c_str());

		/* copy sequence adding $ at beginning and end of sequence (for option aa)*/
		if(addTermini){
			E->S[m][1] = '$';
			for (int pos = 0; pos < E->n-2; pos++) {
				char c = it->seq[pos];
				if (islower(c))	c = static_cast<char>(toupper(c));
				E->S[m][pos + 2] = AlphaCode((int) c, A);
			}
			E->S[m][E->n] = '$';
		/* copy sequence starting at position 1 */
		}else{
			//cerr << E->info[m] << "\t";
			for (int pos = 0; pos < E->n; pos++) {
				char c = it->seq[pos];
				if (islower(c))	c = static_cast<char>(toupper(c));
				E->S[m][pos + 1] = AlphaCode((int) c, A);
				//cerr << c;
			}
			//cerr  << endl;
		}
	}
	return E;
}
