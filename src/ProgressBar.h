#ifndef _PROGRESS_BAR_H
#define _PROGRESS_BAR_H

#include <string>

/**
 * Returns a string à la "[======   ] 75.0%"
 * with total length width.
 */
std::string progressBar(double curr, double total, int width);

#endif
