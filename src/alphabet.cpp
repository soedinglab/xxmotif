#include "alphabet.h"
#include <cctype>
#include <string.h>
#include <stdio.h>

a_type	MkAlpha(const char *map_s)
{
	a_type	A = (alphabet_type*)malloc(sizeof(alphabet_type));
	char	c;
	int n;
	//NEW(A,1,alphabet_type);
	for(n=0;(c=map_s[n])!= '\0';n++) { 			/* Get # letters */
		if(!isalpha(c) && c != '$') {
			fprintf(stderr, "Illegal alphabet string: %c\n", c);
			exit(1);
		}
	}
	A->n = n-1;
	//NEW(A->alphabet,A->n+2,char);				/* ALPHABET */
	A->alphabet = (char*)malloc((A->n+2)*sizeof(char));
	strncpy(A->alphabet,map_s,nAlpha(A)+1);	
	//NEW(A->code2let,(strlen(map_s)+1),char);	/* CODE2LETTER */
	A->code2let = (char*) malloc((strlen(map_s)+1)*sizeof(char));
	strcpy(A->code2let,map_s); 
	
	//NEW(A->let2code,127, int);	/* LETTER2CODE */
	A->let2code = (uint8_t*)malloc(127*sizeof(uint8_t));
	for(uint8_t i=0;i<127;i++)
		A->let2code[i]= 0;	/* =error */
	for(uint8_t i=0;map_s[i]!=0;i++) {
		c = map_s[i]; 
		A->let2code[(int)c] = i; 
		if(isupper(c)) { 
			c = static_cast<char>(tolower(c));
			A->let2code[(int)c] = i; 
		}
		else if(islower(c)) { 
			c = static_cast<char>(toupper(c));
			A->let2code[(int)c] = i; 
		}
	}
	A->C = NULL;
	return (A);
}

char getIUPAC_char(char* string, a_type A){	
	char res = 'N';
	if(nAlpha(A) != 4)	return res;
	if(!strcmp(string, "A")) res = 'A';
	else if(!strcmp(string, "C")) res = 'C';
	else if(!strcmp(string, "G")) res = 'G';
	else if(!strcmp(string, "T")) res = 'T';
	else if(!strcmp(string, "AG")) res = 'R';
	else if(strcmp(string, "CT") == 0) res = 'Y';
	else if(strcmp(string, "CG") == 0) res = 'S';
	else if(strcmp(string, "AT") == 0) res = 'W';
	else if(strcmp(string, "GT") == 0) res = 'K';
	else if(strcmp(string, "AC") == 0) res = 'M';
	else if(strcmp(string, "CGT") == 0) res = 'B';
	else if(strcmp(string, "AGT") == 0) res = 'D';
	else if(strcmp(string, "ACT") == 0) res = 'H';
	else if(strcmp(string, "ACG") == 0) res = 'V';	
	return res; 	
}

void	NilAlpha(a_type A)
{
	if(A != NULL) {
	   free(A->alphabet);
	   free(A->code2let);
	   free(A->let2code);
	   if(A->C != NULL) free(A->C);
	   free(A);
	}
}
