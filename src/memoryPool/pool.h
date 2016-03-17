#ifndef POOL_H_
#define POOL_H_

#include <iostream>
#include <stdlib.h>

class Pool {
public:
	//Pool(unsigned sz) : esize(sz < sizeof(Link*) ? sizeof(Link*) : sz){//; // n is the size of elements
	//	head = 0;
	//	chunks = 0;
	//}
	Pool(long unsigned sz);
	~Pool();
	void *alloc();
	void free(void *b); // put an element back into the pool
	void reset_Pool();
private:
	struct Link {
		Link *next;
	};
	struct Chunk {
		enum {size = 8 * 1024 - 16 };
		//enum {size = 8 * 1048576 - 16 };
		Chunk *next;
		char mem[size];
	};

	Chunk *chunks;
	const long unsigned int esize;
	Link *head;

	Pool(Pool &); // copy protection
	void operator=(Pool &); // copy protection
	void grow();
};

inline Pool::Pool(long unsigned sz) : esize(sz < sizeof(Link*) ? sizeof(Link*) : sz) {
	head = 0;
	chunks = 0;
}

inline Pool::~Pool() // free all chunks
{
	Chunk *n = chunks;
	while (n) {
		Chunk* p = n;
		n = n->next;
		delete p;
	}
}

inline void Pool::reset_Pool(){
	Chunk *n = chunks;
	int nb = 0;
	while (n) {
		Chunk* p = n;
		n = n->next;
		delete p;
		nb++;
	}
//	if(nb > 0){
//		std::cerr << std::endl << nb << " of chunks deleted" << std::endl;
//		std::cerr << "chunksize: " << Chunk::size << "\t" << sizeof(Chunk) << "\nelement size: " << esize <<
//			"\nnb of elements in chunk: " << Chunk::size / esize <<
//			"\n => " << (Chunk::size / 1024.0) * (nb / 1024.0) << " MB freed" << std::endl;
//	}

	head = 0;
	chunks = 0;

	//grow();
}

inline void Pool::grow() // allocate new ‘chunk,’ organize it as a linked list of elements of size ’esize’
{
	Chunk *n = new Chunk;
	n->next = chunks;
	chunks = n;
	const long int nelem = Chunk::size / esize;
	char *start = n->mem;
	char *last = &start[(nelem - 1) * esize];
	for (char *p = start; p < last; p += esize) // assume sizeof(Link)<=esize
		reinterpret_cast<Link*> (p)->next = reinterpret_cast<Link*> (p + esize);
	reinterpret_cast<Link*> (last)-> next = 0;
	head = reinterpret_cast<Link*> (start);
}

inline void* Pool::alloc() {
	if (head == 0)
		grow();
	Link *p = head; // return first element
	head = p->next;
	return p;
}

inline void Pool::free(void *b) {
	Link *p = static_cast<Link*> (b);
	p->next = head; // put b back as first element
	head = p;
}

#endif /* POOL_H_ */
