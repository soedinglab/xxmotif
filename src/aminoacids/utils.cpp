#include "../utils.h" // for log_binomial
#include "utils.h"
#include <cstring>

std::string getBasename(const char* in) {
	char* f = strdup(in);
	char* lastpoint = strrchr(f, '.');
	std::string ret;
	if (lastpoint != NULL) {
		char* n = new char[strlen(f)+1];
		for (size_t i = 0; i <= strlen(f); ++i) {
			n[i] = '\0';
		}
		for (int i = 0; f + i != lastpoint; ++i) {
			n[i] = f[i];
		}
		ret = std::string(n);
		delete []n;
	} else {
		ret = std::string(f);
	}
	free(f);
	return ret;
}
