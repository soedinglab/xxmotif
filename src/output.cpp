#include "output.h"
#include "em/hoUtils.h"

void Output::printTopN(const MotifContainer& motifs, int N){
	int counter = 1;
	for(list<Motif*>::const_iterator it = motifs._startModels.begin();
		it != motifs._startModels.end() && counter <= N; counter++, it++){
		std::cout << counter << ") " << **it;
	}
}

void Output::printMergedMotifs(const MotifContainer& motifs, int N){
	printf("\n\nTOP %d Motifs:\n", N);
	int offset = Global::posSet->max_leng-Global::downstream;

	std::iostream::fmtflags flags = cout.flags();
	std::streamsize prec = cout.precision();

	int printCount = 1;
	for(list<Motif*>::const_iterator it = motifs._startModels.begin();
			it != motifs._startModels.end() && printCount <= N;
			printCount++, it++){
		Motif& m = **it;
//		cerr << m;

		std::string iupac = m.getIUPACString();
		//iupac.erase( remove( iupac.begin(), iupac.end(), ' ' ), iupac.end() );

		const region& r = m.getEnrichment();
		cout << " " << std::left << std::setw(3) << printCount << ": " << std::left << std::setw(34) << iupac << " ";
		cout << "sites: " << std::setw(3) << m.getTotalSites() << std::scientific << std::setprecision(2);
		if(Global::usePositionalProbs){
			if(r.set == 0){
				cout <<	std::right << "  max: " << std::setw(3) << "-" << "\tRegion: -/-";
			}else{
				cout <<	std::right << "  max: " << std::setw(3) << r.max+m.getFirstMotifColumn()-1-offset << "\tRegion: " <<
						r.startRegion+m.getFirstMotifColumn()-1-offset << "/" << r.endRegion+m.getFirstMotifColumn()-1-offset;
			}
		}
		cout <<	"\t\tlog(E-Value): " << m.getPval() << "\tE-Value: ";
		double pVal_log10 = m.getPval() / LogTable::LOG_10;
		if (pow(10, pVal_log10) == 0){
			cout << strprintf("1e%d\t\t ", (int)pVal_log10);
		}else{
			cout << strprintf("%.2e\t\t ", pow(10, pVal_log10) );
		}
		cout << "\tPos Pval: " << pow(10, m.getPosPval()) << endl;
	}
	cout.flags(flags);
	cout << std::setprecision(static_cast<int>(prec));

}

void Output::printFinalMotifs(const MotifContainer& motifs, int N){
	int printCount = 1;
	for(list<Motif*>::const_iterator it = motifs._startModels.begin();
			it != motifs._startModels.end() && printCount <= N;
			it++, printCount++){

		/* output model on standard output */
		printf("\nModel Number %d from %d:\n\n", printCount, motifs.getMotifNb());

		cout << **it;
	}
}


void Output::printOutput(const MotifContainer& motifs, int N) {

	if( motifs.getMotifNb() == 0 ){
		FILE *filepointer;
		char* filename = ( char* )calloc( 1024, sizeof( char ) );
		sprintf( filename, "%s/%s.log", Global::outputDirectory, baseFileName(
				 Global::posSet->name ) );
		filepointer = fopen( filename, "w" );
		free( filename );
		fprintf( filepointer, "No motifs pass filter\n" );
		fclose( filepointer );
		return;
	}

	printFinalMotifs(motifs, N);

	printMotifDistribution(motifs);
	writeSummaryFile(motifs);
	printPvalFile(motifs);
	writeMemeFile(motifs);
	printSequenceFile();
	printScoreDistribution(motifs);
	printPosSetFreq();	


	/* write benchmark data */
	if(Global::benchmarkFolder != NULL) printBenchmark(motifs);
	if(Global::pwmFolder != NULL) print_pwm_folder(motifs);

	/* write cs stuff if necessary */
	if(Global::csbest > 0) print_cs_output(motifs);

	/* write R code into File */
	writeRFile(motifs);
}

void Output::printPosSetFreq(){
	FILE *fptr;
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/tmp/posSetFreq.dat", Global::outputDirectory);
	fptr = fopen(filename, "w");
	free(filename);
	for(int i=1; i<=nAlpha(Global::A); i++){
		fprintf(fptr, "%c\t%.4f\n", AlphaChar(i, Global::A), Global::posBg[i]);
	}
	fclose(fptr);
}


void Output::printBenchmark(const MotifContainer& motifs) {
	char* filename = (char*) calloc(1024,sizeof(char));
	sprintf(filename, "%s_pred.txt", Global::benchmarkFolder);
	FILE* fptr = fopen (filename, "w");
	free(filename);

	Motif& m = **(motifs._startModels.begin());

	for (StartPosContainer::const_iterator it_startPos= m.getStartPosList().begin();
			it_startPos != m.getStartPosList().end();	it_startPos++){
		int seq = it_startPos->seq;
		char* info = (char*)calloc(1024, sizeof(char));
		int i=0;
		while(Global::posSet->entity[seq]->info[0][i] != '\0'){
			info[i]=Global::posSet->entity[seq]->info[0][i];
			i++;
		}
		info[i]='\0';
		fprintf(fptr, "%s, %d, %d, ", info, it_startPos->pos + m.getFirstMotifColumn()-1, it_startPos->pos + m.getLastMotifColumn()-1);
		free(info);
		for(i=it_startPos->pos + m.getFirstMotifColumn()-1; i <= it_startPos->pos + m.getLastMotifColumn()-1; i++){
			fprintf(fptr, "%c", AlphaChar(Global::posSet->entity[seq]->S[0][i] , Global::A));
		}
		fprintf(fptr, "\n");
	}
	fclose(fptr);
}

void Output::print_pwm_folder(const MotifContainer &motifs){
	FILE *fptr = NULL;
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/%s.pwm", Global::pwmFolder, Global::shortFileName);
	fptr = fopen(filename, "w");
	if (fptr == NULL) {
		fprintf(stderr, "ERROR opening file %s for writing PWMs\n", filename);
		return;
	}
	free(filename);
	int motifNb = 1;
	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		double **pwm = (*it_motif)->getPWM();
		fprintf(fptr, "Motif %d: %s", motifNb, (*it_motif)->getIUPACString().c_str());
		double pVal_log10 = (*it_motif)->getPval() / LogTable::LOG_10;
		if (pow(10, pVal_log10) == 0){
			fprintf(fptr, "E-Value: 1e%d\t\t ", (int)pVal_log10);
		}else{
			fprintf(fptr, "E-Value: %.2e\t\t ", pow(10, pVal_log10) );
		}
		fprintf(fptr, "\n");

		for(int i=1;i<=nAlpha(Global::A);i++){
			//fprintf(fptr, "#%c\t", AlphaChar(i,Global::A));
			for(int j = (*it_motif)->getFirstMotifColumn(); j <= (*it_motif)->getLastMotifColumn(); j++){
				fprintf(fptr, "%.5f\t", exp(pwm[j][i]));
			}
			fprintf(fptr, "\n");
		}
		fprintf(fptr, "\n");
	}
	fclose(fptr);
}

void Output::print_cs_output(const MotifContainer& motifs) {
	if(Global::csprofiles.length() == 0){
		cerr << "!!! no output file for csProfiles given !!! \n";
		exit(-1);
	}

	FILE* fptr = fopen (Global::csprofiles.c_str(), "w");
	fprintf(stderr, "print %s\n", Global::csprofiles.c_str());
	double log2 = log(2);
	int motifNb = std::min((int)motifs._startModels.size(), Global::csbest);
	int contextNb = 0;
	list<Motif*>::const_iterator it_motif = motifs._startModels.begin();
	for(int i=0; i<motifNb; i++, it_motif++){
		for(int central_col = (*it_motif)->getFirstMotifColumn(); central_col <= (*it_motif)->getLastMotifColumn(); central_col++){
			if(central_col - Global::cswlen/2 < 1 || central_col + Global::cswlen/2 > PWM_LENGTH) continue;
			contextNb++;
		}
	}
	fprintf(fptr, "ContextLibrary\n");
	fprintf(fptr, "SIZE\t%d\n", contextNb);
	fprintf(fptr, "LENG\t%d\n", Global::cswlen);

	it_motif = motifs._startModels.begin();
	for(int i=0; i<motifNb; i++, it_motif++){
		double **pwm = (*it_motif)->getPWM();
		for(int central_col = (*it_motif)->getFirstMotifColumn(); central_col <= (*it_motif)->getLastMotifColumn(); central_col++){
			if(central_col - Global::cswlen/2 < 1 || central_col + Global::cswlen/2 > PWM_LENGTH) continue;

			fprintf(fptr, "ContextProfile\n");
			string s = (*it_motif)->getIUPACString(NULL);
			s.erase(remove_if(s.begin(), s.end(), isspace), s.end());
			fprintf(fptr, "NAME\t%s%d\n", s.c_str(), central_col-(*it_motif)->getFirstMotifColumn()+1);
			fprintf(fptr, "PRIOR\t%.10e\n", 1.0 / contextNb);
			fprintf(fptr, "ISLOG\tF\n");
			fprintf(fptr, "LENG\t%d\n", Global::cswlen);
			fprintf(fptr, "ALPH\t%u\n", nAlpha(Global::A));

			fprintf(fptr, "\t");
			for(int i=0; i<nAlpha(Global::A); i++)	fprintf(fptr, "%c\t", AlphaChar(i+1, Global::A));
			fprintf(fptr, "\n");

			for(int col = central_col - Global::cswlen/2, j=1; col <= central_col + Global::cswlen/2; col++, j++){
				fprintf(fptr, "%d\t", j);
				for(int i=1; i<=nAlpha(Global::A); i++){
					if(pwm[col][i] == -INFINITY){
						fprintf(fptr, "*\t");
					}else{
						int value = static_cast<int>(round(-1000 * pwm[col][i]/log2));
						fprintf(fptr, "%d\t", value);
					}
				}
				fprintf(fptr, "\n");
			}
			fprintf(fptr, "//\n");
		}
	}
	fclose(fptr);
}

void Output::writeRFile(const MotifContainer& motifs) {
	FILE *fptr;
	int length = Global::posSet->max_leng;
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/tmp/plotDistribution.R", Global::outputDirectory);
	fptr = fopen(filename, "w");
	free(filename);
	fprintf(fptr,
	"####################################################\n" \
	"## Adapted from seqLogo package from Bioconductor ## \n" \
	"####################################################\n" \
	"\n" \
	"## get information content profile from PWM\n" \
	"pwm2ic<-function(pwm) {\n" \
	"    npos<-ncol(pwm)\n" \
	"    ic<-numeric(length=npos)\n" \
	"    for (i in 1:npos) {\n" \
	"        ic[i]<-2 + sum(sapply(pwm[, i], function(x) {\n" \
	"            if (x > 0) { x*log2(x) } else { 0 }\n" \
	"        }))\n" \
	"    }    \n" \
	"    ic\n" \
	"}\n" \
	"\n" \
	"######################\n" \
	"##\n" \
	"## plot sequence logo\n" \
	"##\n" \
	"######################\n" \
	"\n" \
	"letterA <- function(x.pos,y.pos,ht,wt,xspace,id=NULL){\n" \
	"	\n" \
	"  x <- c(0, 4, 6,10, 8,5.0,2,0,  3.2,3.6,6.4,6.8,3.2)\n" \
	"  y <- c(0,10,10, 0, 0,7.5,0,0,  3.0,4.0,4.0,3.0,3.0)\n" \
	"  x <- 0.1*x\n" \
	"  y <- 0.1*y\n" \
	"  \n" \
	"  x <- x.pos + wt*x + xspace\n" \
	"  y <- y.pos + ht*y\n" \
	"  \n" \
	"  if (is.null(id)){\n" \
	"    id <- c(rep(1,9),rep(2,4))\n" \
	"  }else{\n" \
	"    id <- c(rep(id,9),rep(id+1,4))\n" \
	"  }\n" \
	"  \n" \
	"  fill <- c(\"#008000\",\"#008000\")\n" \
	"  col <- c(\"#008000\",\"#008000\")\n" \
	"  \n" \
	"  list(x=x,y=y,id=id,fill=fill,col=col)\n" \
	"}\n" \
	"\n" \
	"## T\n" \
	"letterT <- function(x.pos,y.pos,ht,wt,xspace,id=NULL){\n" \
	"	\n" \
	"  x <- c(0.2,9.8,9.8,6,6,4,4,0.2)\n" \
	"  y <- c(10,10,9,9,0,0,9,9)\n" \
	"  x <- 0.1*x\n" \
	"  y <- 0.1*y\n" \
	"  \n" \
	"  x <- x.pos + wt*x + xspace\n" \
	"  y <- y.pos + ht*y\n" \
	"  \n" \
	"  if (is.null(id)){\n" \
	"    id <- rep(1,8)\n" \
	"  }else{\n" \
	"    id <- rep(id,8)\n" \
	"  }\n" \
	"  \n" \
	"  fill <- \"#FF0000\"\n" \
	"  col <-  \"#FF0000\"\n" \
	"  \n" \
	"  list(x=x,y=y,id=id,fill=fill,col=col)\n" \
	"}\n" \
	"\n" \
	"## C\n" \
	"letterC <- function(x.pos,y.pos,ht,wt,xspace,id=NULL){\n" \
	"  angle1 <- seq(0.3+pi/2,pi,length=100)\n" \
	"  angle2 <- seq(pi,1.5*pi,length=100)\n" \
	"  x.l1 <- 0.5 + 0.5*sin(angle1)\n" \
	"  y.l1 <- 0.5 + 0.5*cos(angle1)\n" \
	"  x.l2 <- 0.5 + 0.5*sin(angle2)\n" \
	"  y.l2 <- 0.5 + 0.5*cos(angle2)\n" \
	"  \n" \
	"  x.l <- c(x.l1,x.l2)\n" \
	"  y.l <- c(y.l1,y.l2)\n" \
	"  \n" \
	"  x <- c(x.l,rev(x.l))\n" \
	"  y <- c(y.l,1-rev(y.l))\n" \
	"  \n" \
	"  x.i1 <- 0.5 +0.35*sin(angle1)\n" \
	"  y.i1 <- 0.5 +0.35*cos(angle1)\n" \
	"  x.i1 <- x.i1[y.i1<=max(y.l1)]\n" \
	"  y.i1 <- y.i1[y.i1<=max(y.l1)]\n" \
	"  y.i1[1] <- max(y.l1)\n" \
	"  \n" \
	"  x.i2 <- 0.5 +0.35*sin(angle2)\n" \
	"  y.i2 <- 0.5 +0.35*cos(angle2)\n" \
	"  \n" \
	"  x.i <- c(x.i1,x.i2)\n" \
	"  y.i <- c(y.i1,y.i2)\n" \
	"  \n" \
	"  x1 <- c(x.i,rev(x.i))\n" \
	"  y1 <- c(y.i,1-rev(y.i))\n" \
	"  \n" \
	"  x <- c(x,rev(x1))\n" \
	"  y <- c(y,rev(y1))\n" \
	"  \n" \
	"  x <- x.pos + wt*x + xspace\n" \
	"  y <- y.pos + ht*y\n" \
	"  \n" \
	"  if (is.null(id)){\n" \
	"    id <- rep(1,length(x))\n" \
	"  }else{\n" \
	"    id <- rep(id,length(x))\n" \
	"  }\n" \
	"  \n" \
	"  fill <- \"#0000FF\"\n" \
	"  col <- \"#0000FF\"\n" \
	"  \n" \
	"  list(x=x,y=y,id=id,fill=fill,col=col)\n" \
	"  }\n" \
	"  \n" \
	"  \n" \
	"## G\n" \
	"letterG <- function(x.pos,y.pos,ht,wt,xspace,id=NULL){\n" \
	"  angle1 <- seq(0.3+pi/2,pi,length=100)\n" \
	"  angle2 <- seq(pi,1.5*pi,length=100)\n" \
	"  x.l1 <- 0.5 + 0.5*sin(angle1)\n" \
	"  y.l1 <- 0.5 + 0.5*cos(angle1)\n" \
	"  x.l2 <- 0.5 + 0.5*sin(angle2)\n" \
	"  y.l2 <- 0.5 + 0.5*cos(angle2)\n" \
	"  \n" \
	"  x.l <- c(x.l1,x.l2)\n" \
	"  y.l <- c(y.l1,y.l2)\n" \
	"  \n" \
	"  x <- c(x.l,rev(x.l))\n" \
	"  y <- c(y.l,1-rev(y.l))\n" \
	"  \n" \
	"  x.i1 <- 0.5 +0.35*sin(angle1)\n" \
	"  y.i1 <- 0.5 +0.35*cos(angle1)\n" \
	"  x.i1 <- x.i1[y.i1<=max(y.l1)]\n" \
	"  y.i1 <- y.i1[y.i1<=max(y.l1)]\n" \
	"  y.i1[1] <- max(y.l1)\n" \
	"  \n" \
	"  x.i2 <- 0.5 +0.35*sin(angle2)\n" \
	"  y.i2 <- 0.5 +0.35*cos(angle2)\n" \
	"  \n" \
	"  x.i <- c(x.i1,x.i2)\n" \
	"  y.i <- c(y.i1,y.i2)\n" \
	"  \n" \
	"  x1 <- c(x.i,rev(x.i))\n" \
	"  y1 <- c(y.i,1-rev(y.i))\n" \
	"  \n" \
	"  x <- c(x,rev(x1))\n" \
	"  y <- c(y,rev(y1))\n" \
	"  \n" \
	"  h1 <- max(y.l1)\n" \
	"  r1 <- max(x.l1)\n" \
	"  \n" \
	"  h1 <- 0.4\n" \
	"  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)\n" \
	"  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)\n" \
	"  \n" \
	"  \n" \
	"  \n" \
	"  if (is.null(id)){\n" \
	"    id <- c(rep(1,length(x)),rep(2,length(x.add)))\n" \
	"  }else{\n" \
	"    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))\n" \
	"  }\n" \
	"  \n" \
	"  x <- c(rev(x),x.add)\n" \
	"  y <- c(rev(y),y.add)\n" \
	"  \n" \
	"  x <- x.pos + wt*x + xspace\n" \
	"  y <- y.pos + ht*y\n" \
	"  \n" \
	"  \n" \
	"  fill <- c(\"#FFA500\",\"#FFA500\")\n" \
	"  col  <- c(\"#FFA500\",\"#FFA500\")\n" \
	"  \n" \
	"  list(x=x,y=y,id=id,fill=fill,col=col)\n" \
	"  \n" \
	"}\n" \
	"\n" \
	"\n" \
	"\n" \
	"addLetter <- function(letters,which,x.pos,y.pos,ht,wt,xspace){\n" \
	"	\n" \
	"  if (which == \"A\"){\n" \
	"    letter <- letterA(x.pos,y.pos,ht,wt,xspace)\n" \
	"  }else if (which == \"C\"){\n" \
	"    letter <- letterC(x.pos,y.pos,ht,wt,xspace)\n" \
	"  }else if (which == \"G\"){\n" \
	"    letter <- letterG(x.pos,y.pos,ht,wt,xspace)\n" \
	"  }else if (which == \"T\"){\n" \
	"    letter <- letterT(x.pos,y.pos,ht,wt,xspace)\n" \
	"  }else{\n" \
	"    stop(\"which must be one of A,C,G,T\")\n" \
	"  }\n" \
	"\n" \
	"  letters$x <- c(letters$x,letter$x)\n" \
	"  letters$y <- c(letters$y,letter$y)\n" \
	"  \n" \
	"  lastID <- ifelse(is.null(letters$id),0,max(letters$id))\n" \
	"  letters$id <- c(letters$id,lastID+letter$id)\n" \
	"  letters$fill <- c(letters$fill,letter$fill)\n" \
	"  letters$col <- c(letters$col,letter$col)\n" \
	"  letters\n" \
	"}\n" \
	"\n" \
	"\n" \
	"\n" \
	"## plot a sequence logo\n" \
	"seqLogo <- function(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15, region=c(27,15,5,15)){\n" \
	"	\n" \
	"  if (class(pwm) == \"pwm\"){\n" \
	"    pwm <- pwm@pwm    \n" \
	"  }else if (class(pwm) == \"data.frame\"){\n" \
	"    pwm <- as.matrix(pwm)\n" \
	"  }else if (class(pwm) != \"matrix\"){\n" \
	"    stop(\"pwm must be of class matrix or data.frame\")\n" \
	"  }\n" \
	"  \n" \
	"  if (any(abs(1 - apply(pwm,2,sum)) > 0.01))\n" \
	"    stop(\"Columns of PWM must add up to 1.0\")\n" \
	"    \n" \
	"    \n" \
	"  chars <- c(\"A\",\"C\",\"G\",\"T\")\n" \
	"  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)\n" \
	"  npos <- ncol(pwm)\n" \
	"  \n" \
	"  \n" \
	"  if (ic.scale){\n" \
	"    ylim <- 2\n" \
	"    ylab <- \"Information content\"\n" \
	"    facs <- pwm2ic(pwm)\n" \
	"  }else{\n" \
	"    ylim <- 1\n" \
	"    ylab <- \"Probability\"\n" \
	"    facs <- rep(1, npos)\n" \
	"  }\n" \
	"  \n" \
	"  xspace <- 0.05\n" \
	"  wt <- 0.95\n" \
	"  x.pos <- 0  \n" \
	"  for (j in 1:npos){\n" \
	"  \n" \
	"    column <- pwm[,j]\n" \
	"    hts <- 0.95*column*facs[j]\n" \
	"    letterOrder <- order(hts)\n" \
	"    \n" \
	"    y.pos <- 0\n" \
	"    for (i in 1:4){\n" \
	"      letter <- chars[letterOrder[i]]\n" \
	"      ht <- hts[letterOrder[i]]\n" \
	"      if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt,xspace)\n" \
	"      y.pos <- y.pos + ht + 0.01\n" \
	"    }\n" \
	"    x.pos <- x.pos + wt\n" \
	"  }\n" \
	"  \n" \
	"  pushViewport(plotViewport(region))\n" \
	"  pushViewport(dataViewport(0:ncol(pwm),0:ylim,name=\"vp1\"))\n" \
	"  \n" \
	"  grid.polygon(x=unit(c(0,ncol(pwm)+1,ncol(pwm)+1,0),\"native\"), y=unit(c(0,0,2,2),\"native\"), \n" \
	"				gp=gpar(fill=\"white\",col=\"white\"))\n" \
	"  grid.polygon(x=unit(letters$x,\"native\"), y=unit(letters$y,\"native\"),\n" \
	"               id=letters$id,\n" \
	"               gp=gpar(fill=letters$fill,col=letters$col,lwd=1))\n" \
	"  if (xaxis){\n" \
	"    grid.xaxis(at=seq(0.5,ncol(pwm)-0.5),label=1:ncol(pwm), gp=gpar(fontsize=xfontsize))\n" \
	"    #grid.text(\"Position\",y=unit(-3,\"lines\"), gp=gpar(fontsize=xfontsize))\n" \
	"  }\n" \
	"  if (yaxis){\n" \
	"     grid.yaxis(gp=gpar(fontsize=yfontsize))\n" \
	"    grid.text(ylab,x=unit(-3,\"lines\"),rot=90, gp=gpar(fontsize=yfontsize))\n" \
	"  }\n" \
	"  popViewport()\n" \
	"  popViewport()\n" \
	"  par(ask=FALSE)\n" \
	"}\n" \
	"\n" \
	"library(grid)\n" \
	"path = commandArgs(TRUE)[1]\n" \
	"data <- read.csv(paste(path, \"/%s.pwm\", sep=\"\"), header=F, sep=\"\\t\", col.names=c(1:30))\n" \
	"\n", Global::shortFileName);

	if(motifs.getMotifNb() != 1){
		fprintf(fptr,
		"motifDistr <- read.csv(paste(path, \"/tmp/motifDistribution.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
		"motifDistr = t(motifDistr) \n" \
		"motifNb = motifDistr[1,]     # extract motif numbers\n" \
		"colnames = motifDistr[2,]     # store colomn names\n" \
		"occ = as.numeric(motifDistr[length(motifDistr[,1])-4,]) # extract occurrence\n" \
		"Pval = as.numeric(motifDistr[length(motifDistr[,1])-3,]) # extract Pvalue\n" \
		"my = (as.numeric(motifDistr[length(motifDistr[,1])-2,]) - %d) # extract my of distribution\n",
			Global::usePositionalProbs ? length-Global::downstream : 0 );
		fprintf(fptr,
		"startRegion = as.numeric(motifDistr[length(motifDistr[,1])-1,]) # extract sigma of distribution\n" \
		"endRegion = as.numeric(motifDistr[length(motifDistr[,1]),]) # extract startRange of distribution\n" \
		"motifDistr = motifDistr[c(3:(length(motifDistr[,1])-5)),]  # extract counts\n" \
		"countsForward <- read.csv(paste(path, \"/tmp/countsForward.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
		"countsBackward <- read.csv(paste(path, \"/tmp/countsBackward.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
		"rounds=ceiling(length(motifNb)/12)\n" \
		"smooth <- 3 \n" \
		"sm <- floor(smooth/2)\n" \
		"Height <- 1200\n" \
		"Width <- Height*(16/9)\n" \
		"rows <- 3\n" \
		"cols <- 4\n" \
		"nb <- Height/15\n" \
		"plotWidth=(nb/cols)*(16/9)+1.1\n" \
		"plotHeight=nb/rows + 1.1\n" \
		"for(offset in c(0:(rounds-1))*rows*cols){\n" \
		"   png(paste(path, \"/%s_\",offset,\".png\",sep=\"\"),width=Width,height=Height)\n",
			Global::shortFileName);
		fprintf(fptr,
		"   colnam = colnames[c((1+offset):(12+offset))]\n" \
		"   layout(matrix(c(1:(rows*cols*2)),rows*2,cols), heights = rep(1,rows*2))\n" \
		"	for (i in (1+offset):(min(length(motifNb), length(colnam)+offset))){\n" \
		"      row=(i-offset-1)%%%%rows \n" \
		"      col=floor((i-offset-1)/rows) \n" \
		"	   par(mai=c(1,1,1,1)/2, mar=c(0.1,2.5,2.5,2.5))\n" \
		"	   a = motifDistr[,i]\n" \
		"	   a = as.numeric(a)\n" \
		"	   b<-barplot(a, cex.main=2, cex.axis=1, main=paste(\"Motif \", motifNb[i], sep=\"\"));\n" \
		"	   bwidth<-(b[2]-b[1])/2\n" \
		"	   lines(c(b[max(1,startRegion[i])]-bwidth,b[min(length(b),endRegion[i])]+bwidth), c(0,0), lwd=10, col=\"red\", type=\"l\")\n");
		if(Global::revcomp){
//			fprintf(fptr, "	lines(c(b[max(1,startRegion[i])]-bwidth,b[min(length(b),endRegion[i])]+bwidth), c(0,0), lwd=10, col=\"red\", type=\"l\")\n");

		}
		fprintf(fptr, "	par(new=T)\n");
		if(Global::usePositionalProbs){
			fprintf(fptr,
			"	   if(Pval[i] < -300)leg.txt=c(paste(\"occ: \",round(occ[i]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: 1e\",round(Pval[i]),sep=\"\"))" \
			"	   else leg.txt=c(paste(\"occ: \",round(occ[i]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: \",format(10**Pval[i], digits=3),sep=\"\"))\n" \
			"	   if(is.finite(my[i])) leg.txt=c(leg.txt, paste(\"max: \",my[i],sep=\"\"))\n");
		}else{
			fprintf(fptr,
			"	   if(Pval[i] < -300)leg.txt=c(paste(\"occ: \",round(occ[i]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: 1e\",round(Pval[i]),sep=\"\"))" \
			"	   else leg.txt=c(paste(\"occ: \",round(occ[i]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: \",format(10**Pval[i], digits=3),sep=\"\"))\n");
		}
		fprintf(fptr,
		"	   legend(\"topleft\", inset=0.0,legend=leg.txt, text.col=\"black\", cex=1.5, box.col=\"black\", box.lwd=1, bg=\"white\", text.width=strwidth(\"babackward backward strand\"), xpd=T)\n" \
		"      pwm <- as.matrix(data[c(2:5)+((i-1)*5),1:min(which(is.na(data[2+(i-1)*5,]))-1)])\n" \
		"      pwm <- matrix(as.double(pwm), nrow=4)\n" \
		"      logoWidth <- 0.3 * sqrt(length(pwm[1,]) / 12)\n" \
		"	   logoCoords <- c(plotHeight*(0.85+rows-1-row), plotWidth*((0.5-logoWidth/2)+col), plotHeight*(0.06+row), plotWidth*((0.5-logoWidth/2)+cols-1-col))\n" \
		"	   seqLogo(pwm, xaxis=F, yaxis=F, xfontsize=10, yfontsize=10, region=logoCoords)\n" \
		"	   par(new=F)\n" \
		"	   par(mar=c(3.5,2.5,0.5,2.5))\n" \
		"	   d <- as.numeric(countsForward[i,])\n" \
		"	   m <- sapply((1+sm):(length(d)-sm), function(a){mean(d[(a-sm):(a+sm)])})\n" \
		"	   plot(c((1+sm):(length(d)-sm))+(%d),m, xlim=c(1,length(b))+(%d),type=\"l\", col=\"red\", lwd=2, cex.axis=1, xlab=\"\")\n",
			Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0, Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0);
		fprintf(fptr,
		"	   d <- as.numeric(countsBackward[i,])\n" \
		"	   m <- sapply((1+sm):(length(d)-sm), function(a){mean(d[(a-sm):(a+sm)])})\n" \
		"	   lines(c((1+sm):(length(d)-sm))+(%d),m, xlim=c(1,length(b))+(%d),type=\"l\", col=\"blue\", lwd=2)\n",
			Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0, Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0);
		fprintf(fptr,
		"	   legend(\"topleft\", c(\"forward strand\", \"backward strand\"), fill=c(\"red\", \"blue\"), inset=0.0, cex=1.5, box.col=\"black\", box.lwd=1, bg=\"white\", xpd=T)\n" \
		"	}\n\n" \
		"   dev.off()\n" \
		"}\n");
	}else{
		fprintf(fptr,
			"Height <- 500\n" \
			"scaleFactor=1.2\n" \
			"Width <- Height*scaleFactor\n" \
			"rows <- 1\n" \
			"cols <- 1\n" \
			"nb <- Height/15\n" \
			"plotWidth=(nb/cols)*scaleFactor+1.1\n" \
			"plotHeight=nb/rows + 1.1\n" \
			"png(paste(path, \"/%s.png\",sep=\"\"),width=Width,height=Height)\n", Global::shortFileName);
		fprintf(fptr,
			"layout(matrix(c(1,2),2,1), heights=c(1,1))\n" \
			"row=0\n" \
			"col=0\n" \
			"par(mai=c(1,1,1,1)/2, mar=c(0.1,2.5,2.5,2.5))\n" \
			"motifDistr <- read.csv(paste(path, \"/tmp/motifDistribution.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
			"motifDistr = t(motifDistr) # sort by pvalue and transpose matrix\n" \
			"motifNb = motifDistr[1,]     # extract motif numbers\n" \
			"colnames = motifDistr[2,]     # store colomn names\n" \
			"occ = as.numeric(motifDistr[length(motifDistr[,1])-4,]) # extract occurrence\n" \
			"Pval = as.numeric(motifDistr[length(motifDistr[,1])-3,]) # extract Pvalue\n" \
			"my = (as.numeric(motifDistr[length(motifDistr[,1])-2,]) - %d) # extract my of distribution\n",
				Global::usePositionalProbs ? length-Global::downstream : 0);
		fprintf(fptr,
			"startRegion = as.numeric(motifDistr[length(motifDistr[,1])-1,]) # extract sigma of distribution\n" \
			"endRegion = as.numeric(motifDistr[length(motifDistr[,1]),]) # extract startRange of distribution\n" \
			"motifDistr = motifDistr[c(3:(length(motifDistr[,1])-5)),]  # extract counts\n" \
			"colnam = colnames[1]\n" \
			"smooth <- 3 \n" \
			"sm <- floor(smooth/2)\n" \
			"a = as.numeric(motifDistr)\n" \
			"b <- barplot(a)\n" \
			"bwidth<-(b[2]-b[1])/2\n" \
			"lines(c(b[max(1,startRegion[1])]-bwidth,b[min(length(b),endRegion[1])]+bwidth), c(0,0), lwd=10, col=\"red\", type=\"l\")\n" \
			"par(new=T)\n");
		if(Global::usePositionalProbs){
			fprintf(fptr, "if(Pval[1] < -300)leg.txt=c(paste(\"occ: \",round(occ[1]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: 1e\",round(Pval[1]),sep=\"\")) ");
			fprintf(fptr, "else leg.txt=c(paste(\"occ: \",round(occ[1]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: \",format(10**Pval[1], digits=3),sep=\"\"))\n");
			fprintf(fptr, "if(is.finite(my[1])) leg.txt=c(leg.txt, paste(\"max: \",my[1],sep=\"\"))\n");
		}else{
			fprintf(fptr, "if(Pval[1] < -300)leg.txt=c(paste(\"occ: \",round(occ[1]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: 1e\",round(Pval[1]),sep=\"\")) ");
			fprintf(fptr, "else leg.txt=c(paste(\"occ: \",round(occ[1]*100)/100,\"%%\",sep=\"\"),paste(\"E-Val: \",format(10**Pval[1], digits=3),sep=\"\"))\n");
		}
		fprintf(fptr,
			"legend(\"topleft\", inset=0.0,legend=leg.txt, text.col=\"black\", cex=1.0, box.col=\"black\", box.lwd=1, bg=\"white\", text.width=strwidth(\"i backward strand\"), xpd=T)\n" \
			"pwm <- as.matrix(data[c(2:5),1:min(which(is.na(data[2,]))-1)])\n" \
			"pwm <- matrix(as.double(pwm), nrow=4)\n" \
			"logoWidth <- 0.3 * sqrt(length(pwm[1,]) / 12)\n" \
			"logoCoords <- c(plotHeight*(0.85+rows-1-row), plotWidth*((0.5-logoWidth/2)+col), plotHeight*(0.06+row), plotWidth*((0.5-logoWidth/2)+cols-1-col))\n" \
			"seqLogo(pwm, xaxis=F, yaxis=F, xfontsize=10, yfontsize=10, region=logoCoords)\n" \
			"par(new=F)\n" \
			"par(mar=c(2.5,2.5,0.5,2.5))\n" \
			"countsForward <- read.csv(paste(path, \"/tmp/countsForward.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
			"countsBackward <- read.csv(paste(path, \"/tmp/countsBackward.dat\", sep=\"\"), sep=\"\\t\", header=FALSE)\n" \
			"d <- as.numeric(countsForward[1,])\n" \
			"m <- sapply((1+sm):(length(d)-sm), function(a){mean(d[(a-sm):(a+sm)])})\n" \
			"plot(c((1+sm):(length(d)-sm))+(%d),m, xlim=c(1,length(b))+(%d),type=\"l\", col=\"red\", lwd=2, cex.axis=1, xlab=\"\")\n",
				Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0, Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0);
		fprintf(fptr,
			"d <- as.numeric(countsBackward[1,])\n" \
			"m <- sapply((1+sm):(length(d)-sm), function(a){mean(d[(a-sm):(a+sm)])})\n" \
			"lines(c((1+sm):(length(d)-sm))+(%d),m, xlim=c(1,length(b))+(%d),type=\"l\", col=\"blue\", lwd=2)\n",
				Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0, Global::usePositionalProbs ? (int)(-length+Global::downstream) : 0);
		fprintf(fptr,
			"legend(\"topleft\", c(\"forward strand\", \"backward strand\"), fill=c(\"red\", \"blue\"), inset=0.0, cex=1.0, box.col=\"black\", box.lwd=1, bg=\"white\", xpd=T)\n");

	}
	fclose(fptr);
}

/* TODO: comment (author: siebert)*/
void Output::writeBestKmerSites( MotifContainer& motifs, int minSites, char* baseFileName ){

	Motif* motif;

	int pos, lastPos;
	uint8_t* sequence;

	list<Motif*>::const_iterator iter;
	StartPosContainer::iterator sites;

	std::stringstream str;
	str.str( "" );
	str << Global::outputDirectory << "/" << baseFileName << ".blocks";
	FILE* fptr = fopen( str.str().c_str(), "w" );

	for( iter=motifs.getMotifs().begin(); iter != motifs.getMotifs().end(); iter++ ){

		motif = *iter;

		if( motif->getTotalSites() >= minSites ){

			for( sites=motif->getStartPosList().begin(); sites != motif->getStartPosList().end(); sites++ ){

				sequence = Global::posSet->entity[ sites->seq ]->S[0];
				pos = sites->pos + motif->getFirstMotifColumn() - 1;
				lastPos = pos + motif->getMotifLength() - 1;

				for( ; pos <= lastPos; ++pos )
					fprintf( fptr, "%c", AlphaChar( sequence[pos], Global::A) );
				fprintf( fptr, "\n" );
			}

			printf( "%s\t%d/%d (%.2f%%)\n", motif->getIUPACString().c_str(), motif->getTotalSites(), motif->getPosSetSize(), motif->getTotalSites() * 100.0 / motif->getPosSetSize() );

			break;
		}
	}

	fclose( fptr );
}

void Output::writeBlocksFile(const MotifContainer& motifs){
	FILE *fptr = NULL;
	int motifNb = 0;

	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		char* filename = (char*)calloc(1024, sizeof(char));
		sprintf(filename, "%s/tmp/MotifBlocks%d.blocks", Global::outputDirectory, motifNb);
		fptr = fopen (filename, "w");
		free(filename);

		for (StartPosContainer::const_iterator it_startPos= (*it_motif)->getStartPosList().begin();
				it_startPos != (*it_motif)->getStartPosList().end(); it_startPos++){
			//fprintf(stderr, "%4d/%4d\t", it_startPos->seq, it_startPos->pos + (*it_motif)->getFirstMotifColumn() - 1 - 20);
			for(int l=0; l< 1; l++){ // only for first specie
				for(int k=it_startPos->pos+(*it_motif)->getFirstMotifColumn()-1; k <= it_startPos->pos + (*it_motif)->getLastMotifColumn()-1; k++){
				//for(int k=it_startPos->pos+(*it_motif)->getFirstMotifColumn()-1-20; k <= it_startPos->pos + (*it_motif)->getLastMotifColumn()-1+20; k++){
					if(k < 1 || k > Global::posSet->entity[it_startPos->seq]->n) continue;
					char outputChar = AlphaChar(Global::posSet->entity[it_startPos->seq]->S[l][k], Global::A);
					if(outputChar == '$') outputChar = 'Z';
					fprintf(fptr, "%c", outputChar);
				}
				fprintf(fptr, "\n");
			}
		}
		fclose(fptr);
	}
}

void Output::printMotifDistribution(const MotifContainer& motifs){
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/tmp/motifDistribution.dat", Global::outputDirectory);
	FILE *f_motifDistribution = fopen(filename, "w");
	free(filename);

	int motifNb = 0;
	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		Motif& m = **it_motif;
		std::string title = m.getIUPACString().c_str();

		int* countStartPos = (int*)calloc(Global::posSet->max_leng+1, sizeof(int));
		StartPosContainer& startPosList = m.getStartPosList();
		int totalSites = m.getTotalSites();

		for(StartPosContainer::iterator it = startPosList.begin(); it != startPosList.end(); it++){
			countStartPos[it->pos + m.getFirstMotifColumn() - 1]++;
		}

		fprintf(f_motifDistribution, "%d\t%s\t", motifNb+1, title.c_str());
		for (int j = 1; j <= Global::posSet->max_leng; j++) {
			fprintf(f_motifDistribution, "%d\t", countStartPos[j]);
		}
		free(countStartPos);

//		fprintf(f_motifDistribution, "%f\t%f\t%d\t%d\t%d\n",
//						totalSites * 100.0 / m.getPosSetSize(),
//						m.getPval()/LogTable::LOG_10,
//						m.getEnrichment().max + m.getFirstMotifColumn()-1,
//						m.getEnrichment().startRegion + m.getFirstMotifColumn()-1,
//						m.getEnrichment().endRegion + m.getFirstMotifColumn()-1);

		fprintf(f_motifDistribution, "%f\t%f\t",
				totalSites * 100.0 / m.getPosSetSize(),
				m.getPval()/LogTable::LOG_10);

		if(Global::usePositionalProbs){
			if(m.getEnrichment().set == 0){
				fprintf(f_motifDistribution, "NA\tNA\tNA\n");
			}else{
				fprintf(f_motifDistribution, "%d\t%d\t%d\n",
						m.getEnrichment().max + m.getFirstMotifColumn()-1,
						m.getEnrichment().startRegion + m.getFirstMotifColumn()-1,
						m.getEnrichment().endRegion + m.getFirstMotifColumn()-1);
			}
		}else{
			fprintf(f_motifDistribution, "NA\tNA\tNA\n");
		}
	}

	fclose(f_motifDistribution);
}

void Output::writeSummaryFile(const MotifContainer& motifs){
	int offset = (Global::posSet->max_leng - Global::downstream);

	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/%s_MotifFile.txt", Global::outputDirectory, Global::shortFileName);
	FILE *f_MotifFile = fopen(filename, "w");
	free(filename);

	list<Motif*>::const_iterator it_motif = motifs._startModels.begin();

	fprintf(f_MotifFile, "============= RESULTS MOTIF SEARCH =============\n\n");
	if (Global::startMotif != NULL)	fprintf(f_MotifFile, "start motif: %s\n", Global::startMotif);

	char* IUPAC_long = (char*)calloc(256, sizeof(char));
	for(int motifNb = 0; it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		Motif& m = **it_motif;
		fprintf(f_MotifFile, "\n\t----Motif Number %d: ----\n => ", motifNb+1);

		string IUPAC = m.getIUPACString(IUPAC_long);
		fprintf(f_MotifFile, "%s\n\t", IUPAC_long);

		double pVal_log10 = m.getPval() / LogTable::LOG_10;
		if(pow(10, pVal_log10) == 0){
			fprintf(f_MotifFile, "E-Value: 1e%d\t", (int) pVal_log10 );
		}else{
			fprintf(f_MotifFile, "E-Value: %.2e,\t", pow(10, pVal_log10));
		}
		if(Global::multipleOccurrence){
			fprintf(f_MotifFile, "%d sites in %d sequences\n", m.getTotalSites(), m.getPosSetSize() );
		}else{
			fprintf(f_MotifFile, "occ: %.2f %% (%d of %d sequences)\n", m.getTotalSites() * 100.0 / m.getPosSetSize(), m.getTotalSites(), m.getPosSetSize() );
		}
		fprintf(f_MotifFile, "\tIUPAC: %s\n", IUPAC.c_str());
		if(Global::usePositionalProbs){
			if(m.getEnrichment().set == 0){
				fprintf(f_MotifFile, "\tmode: -, startRegion: -\t endRegion: -\n");
			}else{
				fprintf(f_MotifFile, "\tmode: %d, startRegion: %d\t endRegion: %d\n",
						m.getEnrichment().max + m.getFirstMotifColumn()-1 - offset,
						m.getEnrichment().startRegion + m.getFirstMotifColumn()-1 - offset,
						m.getEnrichment().endRegion + m.getFirstMotifColumn()-1 - offset);
			}
		}
	}
	free(IUPAC_long);
	fclose(f_MotifFile);
}

void Output::printStartPosFile(const MotifContainer& motifs){
	int offset = (Global::posSet->max_leng - Global::downstream);

	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/%s_StartPos.txt", Global::outputDirectory,
			Global::shortFileName);
	FILE *f_startPos = fopen(filename, "w");
	free(filename);

	int motifNb = 0;
	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		Motif& m = **it_motif;
		fprintf(f_startPos, "%d\t%s\t", motifNb, m.getIUPACString().c_str());

		for (StartPosContainer::const_iterator it_startPos= m.getStartPosList().begin();
			it_startPos != m.getStartPosList().end(); it_startPos++){
			//fprintf(f_startPos, "%s,%d,%d;", Global::posSet->entity[it_startPos->seq]->info[0], it_startPos->seq, it_startPos->pos + m.getFirstMotifColumn() -1 - offset);
			fprintf(f_startPos, "%d,%d;", it_startPos->seq, it_startPos->pos + m.getFirstMotifColumn() -1 - offset);

		}
		fprintf(f_startPos, "\n");
	}

	fclose(f_startPos);
}

void Output::printPvalFile(const MotifContainer& motifs){

	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/%s_Pvals.txt", Global::outputDirectory,
			Global::shortFileName);
	FILE *f_startPos = fopen(filename, "w");
	free(filename);

	int motifNb = 1;
	int colLength = 142;
	if(Global::posSet->max_MultSeq > 1) colLength += 16;
	for(int i=0; i<colLength; i++)	fprintf(f_startPos, "#");
	fprintf(f_startPos, "\n###%22s\t%20s\t%20s\t%10s\t%10s\t%10s\t%10s","upstream site", "site", "downstream site", "seq nb","start pos", "strand", "site Pval");
	if(Global::posSet->max_MultSeq > 1)fprintf(f_startPos, "\t%10s", "cons Pval");
	fprintf(f_startPos, " ###\n");
	for(int i=0; i<colLength; i++)	fprintf(f_startPos, "#");
	fprintf(f_startPos, "\n\n");

	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		Motif* m = *it_motif;

		double combined_pValCons;

		list<double> pValCons;	list<double> pValSite;
		char kmerLeft[30];	char kmer[20];	char kmerRight[30];

		StartPosUpdater::getInstance().getProbabilitiesPerSite(pValSite, pValCons, combined_pValCons, m);

		fprintf(f_startPos, "Motif %d: %s\t", motifNb, m->getIUPACString().c_str());
		double pVal_log10 = (*it_motif)->getPval() / LogTable::LOG_10;
		if (pow(10, pVal_log10) == 0){
			fprintf(f_startPos, "E-Value: 1e%d", (int)pVal_log10);
		}else{
			fprintf(f_startPos, "E-Value: %.2e", pow(10, pVal_log10) );
		}
		if(Global::posSet->max_MultSeq > 1)	fprintf(f_startPos, "\t(cons P-value: %.2e)", exp(combined_pValCons));
		fprintf(f_startPos, "\n");

		StartPosContainer::const_iterator it_startPos= m->getStartPosList().begin();
		list<double>::const_iterator it_cons = pValCons.begin();
		list<double>::const_iterator it_site = pValSite.begin();
		for(; it_cons != pValCons.end(); it_cons++, it_site++, it_startPos++){
			int seq = it_startPos->seq;
			int pos = it_startPos->pos + m->getFirstMotifColumn() - 1;

			char strand = '+';
			int seqLength = Global::posSet->entity[seq]->n;
			if(Global::revcomp) seqLength = (seqLength - 1)/2;

			int j=0;
			int start = std::max(1,pos-20);
			if(Global::revcomp && pos > seqLength) start = std::max(seqLength+2, pos-20);
			for(int i=start; i<pos; i++, j++){
				kmerLeft[j] = AlphaChar(Global::posSet->entity[seq]->S[0][i], Global::A);
			}
			kmerLeft[j] = '\0';

			for(int i=0; i<m->getMotifLength(); i++){
				kmer[i] =  AlphaChar(Global::posSet->entity[seq]->S[0][pos+i], Global::A);
			}
			kmer[m->getMotifLength()] = '\0';

			j=0;
			int end = std::min(Global::posSet->entity[seq]->n+1, pos+m->getMotifLength() + 20);
			if(Global::revcomp && pos < seqLength) end = std::min(seqLength+1, pos+m->getMotifLength() + 20);
			for(int i=pos+m->getMotifLength(); i< end; i++, j++){
				kmerRight[j] = AlphaChar(Global::posSet->entity[seq]->S[0][i], Global::A);
			}
			kmerRight[j] = '\0';

			if(Global::revcomp && pos > seqLength){
				strand = '-';
				pos = Global::posSet->entity[seq]->n - (pos + m->getMotifLength() - 1) + 1;
				//for(int i=0; i<m->getMotifLength(); i++){
				//	fprintf(stderr, "%c", AlphaChar(Global::posSet->entity[seq]->S[0][pos+i], Global::A));
				//}
				//fprintf(stderr, "\n");
			}
			fprintf(f_startPos, "%25s\t%20s\t%20s\t%10d\t%10d\t%10c\t%10.4e",
					kmerLeft, kmer, kmerRight, seq, pos, strand, *it_site);
			if(Global::posSet->max_MultSeq > 1) fprintf(f_startPos, "\t%10.4e", *it_cons);
			fprintf(f_startPos, "\n");
		}

		fprintf(f_startPos, "\n");
	}

	fclose(f_startPos);
}

void Output::writeMemeFile(const MotifContainer& motifs) {
  char* filename = (char*)calloc(1024, sizeof(char));
  sprintf(filename, "%s/%s.meme", Global::outputDirectory,
      Global::shortFileName);
  FILE *fmeme = fopen(filename, "w");
  free(filename);

  fprintf(fmeme, "MEME version 4\n\n");
  fprintf(fmeme, "ALPHABET= ");
  for (int i=1; i<=nAlpha(Global::A); ++i) {
    fprintf(fmeme, "%c", AlphaChar(i, Global::A));
  }
  fprintf(fmeme, "\n\n");
  fprintf(fmeme, "strands: %s\n\n", Global::revcomp ? "+ -" : "+");
  fprintf(fmeme, "Background letter frequencies\n");
  for (int i=1; i<=nAlpha(Global::A); ++i) {
    if (i>1) fprintf(fmeme, " ");
    fprintf(fmeme, "%c %.3f", AlphaChar(i, Global::A), Global::negSet == NULL ? Global::posBg[i] : Global::negBg[i]);
  }
  fprintf(fmeme, "\n");

  {
    int motcount = 1;
    for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); ++it_motif, ++motcount){
      const Motif * const m = *it_motif;
      const double pval = exp(m->getPval());
      const int nsites = m->getTotalSites();
      const int length = m->getLastMotifColumn() - m->getFirstMotifColumn() + 1;
      const double * const * const pwm = (*it_motif)->getPWM();
      fprintf(fmeme, "\nMOTIF %d\n", motcount);
      fprintf(fmeme, "letter-probability matrix: alength= %d w= %d nsites= %d E= %.2e\n",
          nAlpha(Global::A), length, nsites, pval);
      for (int i=m->getFirstMotifColumn(); i<=(*it_motif)->getLastMotifColumn(); ++i) {
        for (int j=1; j<=nAlpha(Global::A); ++j) {
          if (j>1) fprintf(fmeme, " ");
          fprintf(fmeme, "%9.6f", exp(pwm[i][j]));
        }
        fprintf(fmeme, "\n");
      }
    }
  }

  fclose(fmeme);
}

void Output::printSequenceFile(){

	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/%s_sequence.txt", Global::outputDirectory,
			Global::shortFileName);
	FILE *f_sequence = fopen(filename, "w");
	free(filename);

	int colLength = 60;
	for(int i=0; i<colLength; i++)	fprintf(f_sequence, "#");
	fprintf(f_sequence, "\n###%7s\t%10s\t%17s\t","seq nb","seq length","seq header");
	fprintf(f_sequence, " ###\n");
	for(int i=0; i<colLength; i++)	fprintf(f_sequence, "#");
	fprintf(f_sequence, "\n\n");

	for(int i=1; i<= Global::posSet->nent; i++){
		int seqLength = Global::posSet->entity[i]->n;
		if(Global::revcomp) seqLength = (seqLength - 1)/2;
		char* header = Global::posSet->entity[i]->info[0];
		fprintf(f_sequence, "%10d\t%10d\t%s\n", i, seqLength, header);
	}
	fclose(f_sequence);
}

void Output::printScoreDistribution(const MotifContainer& motifs){
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/tmp/countsForward.dat", Global::outputDirectory);
	FILE *fptr_forward = fopen(filename, "w");
	sprintf(filename, "%s/tmp/countsBackward.dat", Global::outputDirectory);
	FILE *fptr_backward = fopen(filename, "w");
	free(filename);

	int motifNb = 0;

	for(list<Motif*>::const_iterator it_motif = motifs._startModels.begin(); it_motif != motifs._startModels.end(); it_motif++, motifNb++){
		Motif* m = *it_motif;
		double maxScore = m->getMaxMatrixScore();

		Motif* m2 = m->getPalindrome();
		double maxScore2 = m2->getMaxMatrixScore();

		//fprintf(stderr, "maxScore: %f\n", maxScore);
		for (int pos = 1; pos <= Global::posSet->max_leng - m->getUsedLength() + 1; pos++) {
			double counts_forward = 0;
			double counts_backward= 0;
			double scoreThreshold = 0.5;
			for (int seq = 1; seq <= m->getPosSetSize(); seq++) {
				if (Global::posSet->entity[seq]->n - m->getMotifLength() + 1 < pos) continue;
				double score_forward = m->calculateMatrixScore_logOdds(seq, pos);
				double score_backward = m2->calculateMatrixScore_logOdds(seq,pos);
				if(score_forward > maxScore * scoreThreshold) counts_forward++;
				if(score_backward > maxScore2 * scoreThreshold) counts_backward++;
			}
			fprintf(fptr_forward, "%f\t", counts_forward);
			fprintf(fptr_backward, "%f\t", counts_backward);
		}

		delete m2;

		for (int pos = 1; pos < m->getUsedLength()-1; pos++) {
			fprintf(fptr_forward, "0\t");
			fprintf(fptr_backward, "0\t");
		}
		fprintf(fptr_forward, "0\n");
		fprintf(fptr_backward, "0\n");
	}

	fclose(fptr_forward);
	fclose(fptr_backward);
}

void Output::printDistanceDistribution(const MotifContainer& motifs){
	char* filename = (char*)calloc(1024, sizeof(char));
	sprintf(filename, "%s/tmp/distanceDistribution.dat", Global::outputDirectory);
	FILE *f_dist_distribution = fopen(filename, "w");
	Motif& m = **(motifs._startModels.begin());

	int* countStartPos = (int*)calloc(Global::posSet->max_leng+1, sizeof(int));
	for(StartPosContainer::iterator it = m.getStartPosList().begin(); it != m.getStartPosList().end(); it++){
		countStartPos[it->pos + m.getFirstMotifColumn() - 1]++;
	}
	for (int j = 1; j <= Global::posSet->entity[1]->n; j++) {
		fprintf(f_dist_distribution, "%d\t%d\n", j - (Global::posSet->entity[1]->n - Global::downstream), countStartPos[j]);
	}
	free(countStartPos);

	fclose(f_dist_distribution);
	free(filename);
}

