#ifndef SORTED_SITES
#define SORTED_SITES

#include <stdint.h>

struct sorted_sites_type{
   	float	score;
   	int32_t		startPos;
   	int32_t 	sequence;
};
typedef sorted_sites_type* sorted_sites;

namespace SortedSites{
	/* compare function for sorting sortedSites_type arrays */
	inline int compare_scores (const void *a, const void *b)
	{
	  sorted_sites *a1 = (sorted_sites *)a;
	  sorted_sites *b1 = (sorted_sites *)b;
	  double diff = (*a1)->score - (*b1)->score;
	  if(diff > 0.0) return 1;
	  else if(diff < 0.0) return -1;
	  return 0;
	}

	/* compare function for sorting sortedSites_type arrays */
	inline int compare_startPos (const void *a, const void *b)
	{
	  sorted_sites *a1 = (sorted_sites *)a;
	  sorted_sites *b1 = (sorted_sites *)b;
	  double diff = (*a1)->sequence - (*b1)->sequence;
	  if(diff > 0.0) return 1;
	  else if(diff < 0.0) return -1;
	  else{
		  diff = (*a1)->startPos - (*b1)->startPos;
		  if(diff > 0.0) return 1;
		  else if(diff < 0.0) return -1;
	  }
	  return 0;
	}

	/* compare function for sorting sortedSites_type arrays */
	inline int compare_startPos_palin (const void *a, const void *b)
	{
	  sorted_sites *a1 = (sorted_sites *)a;
	  sorted_sites *b1 = (sorted_sites *)b;
	  double diff = (*a1)->sequence - (*b1)->sequence;
	  if(diff > 0.0) return 1;
	  else if(diff < 0.0) return -1;
	  else{
		  diff = (*a1)->startPos - (*b1)->startPos;
		  if(diff < 0.0) return 1;
		  else if(diff > 0.0) return -1;
	  }
	  return 0;
	}
}

#endif
