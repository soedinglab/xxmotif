#include "ProgressBar.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <string.h>
#include <stdio.h>

std::string progressBar(double curr, double total, int width) {
	assert(curr >= 0 && curr <= total && width >= 10);
	char buffer[256];
	double p = curr / total;
	sprintf(buffer, "[%*s] %5.1f%%%c", width - 9, "", p * 100, '\0');
	memset(buffer + 1, '=', static_cast<size_t>(ceil(p * (width - 9))));
	return std::string(buffer);
}
