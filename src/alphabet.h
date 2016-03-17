#ifndef ALPHA
#define ALPHA
#include <stdlib.h>
#include <stdint.h>

typedef struct {
	int     n;			/* number of LETTERS */
	char    *alphabet;		/* ALPHABET */
	char    *code2let;		/* CODE2LETTER */
	uint8_t     *let2code;		/* LETTER2CODE */
	char    *C;				/* complementary bases */
} alphabet_type;
typedef	alphabet_type *a_type;

/******************************** PUBLIC ***********************************/
a_type	MkAlpha(const char *map_S);	/* define alphabet */
void	NilAlpha(a_type A);				/* make alphabet A undefined */
char 	getIUPAC_char(char* string, a_type A);	/*return IUPAC character */

/**************************** alphabet defines ****************************/
#define nAlpha(A)		((A)->n)
#define	AlphaChar(i,A)		((A)->code2let[(i)])
#define	AlphaCode(c,A)		((A)->let2code[(int)(c)])

#endif
