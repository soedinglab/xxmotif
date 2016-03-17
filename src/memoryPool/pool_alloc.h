#ifndef POOL_ALLOC_H
#define POOL_ALLOC_H

#include "pool.h"
#include <cassert>
#include <cstddef>
#include <iostream>

using std::size_t;
using std::ptrdiff_t;

template <class T> class Pool_alloc {
public:
	typedef T value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;

	static int poolSize;

	static void reset_Pool(){
		mem.reset_Pool();
	}

	pointer address(reference r) const {
		return &r;
	}

	const_pointer address(const_reference r) const {
		return &r;
	}

	Pool_alloc() throw(){
		/* create allocater class */
	}

	template <class U> Pool_alloc(const Pool_alloc<U>&)throw(){
	}

	~Pool_alloc() throw (){
		/* free allocater class */
	}

	pointer allocate(size_type n, const void* = 0){ //space for n Ts
		/* get memory from memory pool */
		if (n == 1){
			return static_cast<T*> (mem.alloc());
		}
		else{
			assert(false);
			return 0;
		}
	}

	void deallocate(pointer p, size_type n){ //deallocate n Ts,don’t destroy
		/* free memory from memory pool */
		if (n == 1) {
			mem.free(p);
			return;
		} else {
			assert(false);
		}
	}
	void construct(pointer p, const T& val) {
		/* initialize *p by val */
		new(p) T(val);
	}

	void destroy(pointer p) {
		/* destroy *p but don’t deallocate */
		p->~T();
	}

	template <class U>
	struct rebind {
		typedef Pool_alloc<U> other;
	};	//in effect:typedef Pool_alloc<U> other

	size_type max_size() const throw(){
		return size_t(-1) / sizeof(T);
	}

private:
	static Pool mem; // pool of elements of size of(T)

};

template <class T> Pool Pool_alloc<T>::mem(sizeof(T));

#endif /* POOL_ALLOC_H */
