#include "StateLib.h"
#include "../Globals.h"
#include "MProGlobal.h"
#include <fstream>
#include <sstream>

const Alphabet StateLib::dummyAlphabet;

void StateLib::computeSigStates(const double thresh) {
	significanceThreshold = thresh;
	sigStates.resize(alphabet.size());
	for (int i = 0; i < static_cast<int>(alphabet.size()); ++i) {
		if (i == alphabet.wildcardIndex()) {
			continue;
		}
		sigStates[i].clear();
		double bestMatchValue = -std::numeric_limits<double>::max();
		int bestMatchIndex = std::numeric_limits<int>::max();
		for (int s = 0; s < static_cast<int>(states.size()); ++s) {
			if (states[s][i] > bestMatchValue) {
				bestMatchValue = states[s][i];
				bestMatchIndex = s;
			}
			if (states[s][i] - Global::negBg_log[i] >= thresh) {
				std::pair<double, int> p(states[s][i], s);
				sigStates[i].push_back(p);
			}
		}
		/* if no state reaches the significance threshold, take
		 * the one with the highest score for the current aa
		 */
		if (sigStates[i].size() == 0) {
			cerr << "WARNING: no significant state for " << alphabet.itoc(i) <<
					", using best match " << states[bestMatchIndex] << endl;
			sigStates[i].push_back(std::pair<double,int>(bestMatchValue,bestMatchIndex));
		}
		sigStates[i].sort(pairCmp<double, int, std::greater<double> > ());
	}
}

void StateLib::sortByEntropy() {
	std::sort(states.begin(), states.end());
	for (int i = 0; i < static_cast<int>(states.size()); ++i) {
		states[i].setIndex(i);
	}
}

void StateLib::writeTikz(const int per_line, const std::string filename) const {
	static const ColumnState *currentState;
	FILE *fout = fopen(filename.c_str(), "w");
	fprintf(fout, "{\n");
	fprintf(fout, "\t\\newif\\ifshowaa\\showaatrue\n");
	struct AArep {
		char chr;
		int index;
		std::string color;
		int chemgroup;
		static int comp(const void *a, const void *b) {
			AArep* a1 = (AArep*) a;
			AArep* a2 = (AArep*) b;
			if (a1->chemgroup < a2->chemgroup) {
				return -1;
			} else if (a1->chemgroup > a2->chemgroup) {
				return 1;
			} else {
				if ((*currentState)[a1->index] > (*currentState)[a2->index]) {
					return -1;
				} else if ((*currentState)[a1->index] < (*currentState)[a2->index]) {
					return 1;
				} else {
					if (a1->chr < a2->chr) {
						return -1;
					} else if (a1->chr > a2->chr) {
						return 1;
					} else {
						return 0;
					}
				}
			}
		}
	};
	AArep *acids = new AArep[alphabet.size()];

	for (int i = 1; i < static_cast<int>(alphabet.size()); ++i) {
		const char chr = alphabet.itoc(i);
		AArep aa;
		aa.index = i;
		aa.chr = chr;
		std::stringstream line;
		switch (chr) {
			case 'A':
				aa.chemgroup = 4;
				break;
			case 'C':
				aa.chemgroup = 3;
				break;
			case 'D':
				aa.chemgroup = 5;
				break;
			case 'E':
				aa.chemgroup = 5;
				break;
			case 'F':
				aa.chemgroup = 2;
				break;
			case 'G':
				aa.chemgroup = 4;
				break;
			case 'H':
				aa.chemgroup = 7;
				break;
			case 'I':
				aa.chemgroup = 1;
				break;
			case 'K':
				aa.chemgroup = 8;
				break;
			case 'L':
				aa.chemgroup = 1;
				break;
			case 'M':
				aa.chemgroup = 1;
				break;
			case 'N':
				aa.chemgroup = 6;
				break;
			case 'P':
				aa.chemgroup = 9;
				break;
			case 'Q':
				aa.chemgroup = 6;
				break;
			case 'R':
				aa.chemgroup = 8;
				break;
			case 'S':
				aa.chemgroup = 4;
				break;
			case 'T':
				aa.chemgroup = 4;
				break;
			case 'V':
				aa.chemgroup = 1;
				break;
			case 'W':
				aa.chemgroup = 2;
				break;
			case 'Y':
				aa.chemgroup = 2;
				break;
			case '$':
				aa.chemgroup = 0;
				break;
		}
		switch (chr) {
			case 'W':
			case 'Y':
			case 'F':
				aa.color = std::string("00C000");
				break;
			case 'C':
				aa.color = std::string("FFFF00");
				break;
			case 'D':
			case 'E':
				aa.color = std::string("6080FF");
				break;
			case 'L':
			case 'I':
			case 'V':
			case 'M':
				aa.color = std::string("02FF02");
				break;
			case 'K':
			case 'R':
				aa.color = std::string("FF0000");
				break;
			case 'Q':
			case 'N':
				aa.color = std::string("E080FF");
				break;
			case 'H':
				aa.color = std::string("FF8000");
				break;
			case 'P':
				aa.color = std::string("A0A0A0");
				break;
			case 'G':
			case 'A':
			case 'S':
			case 'T':
				aa.color = std::string("FFFFFF");
				break;
			case '$':
				aa.color = std::string("FF789E");
				break;
		}
		fprintf(fout, "\t\\definecolor{rescol%c}{HTML}{%s}\n", (chr == '$' ? 'Z' : chr), aa.color.c_str());
		acids[i] = aa;
	}

	const float scale = 0.1f;
	const float height = 40 * scale;
	const float width = 4 * scale;
	const float sep = 3 * scale;

	for (int s = 0; s < static_cast<int>(states.size()); ++s) {
		currentState = &(states[s]);
		std::qsort(&(acids[1]), alphabet.size() - 1, sizeof(AArep), AArep::comp);
		float y_curr = 0;
		fprintf(fout, "\t\\expandafter\\newcommand\\csname %s%d\\endcsname{\\begin{tikzpicture}[thick]\n",
				Global::type_of_states == CS_63 ? "cs" : "bs", s);
		for (int i = 1; i < alphabet.size(); ++i) {
			float h = height * static_cast<float>(exp((*currentState)[acids[i].index]));
			fprintf(fout, "\t\t\\draw[fill=rescol%c] (%f,%f) rectangle (%f,%f);\n", acids[i].chr == '$' ? 'Z'
					: acids[i].chr, 0.0, y_curr, width, y_curr + h);
			if (h >= 4 * scale) {
				fprintf(fout, "\t\t\\ifshowaa\\node[anchor=center] at (%f,%f) {%s};\\fi\n", 0.5 * width, y_curr + 0.5
						* h, acids[i].chr == '$' ? "\\$" : std::string(1, acids[i].chr).c_str());
			}
			y_curr += h;
		}
		fprintf(fout, "\t\\end{tikzpicture}}\n");
	}

	fprintf(fout, "\t\\sffamily\\bfseries\n\t\\begin{tikzpicture}[line width=0mm,inner sep=0mm]\n");

	float y = 0;
	float x = 0;
	for (int s = 0; s < static_cast<int>(states.size()); ++s) {
		if ((s > 0) && (s % per_line == 0)) {
			x = 0;
			y -= height + sep;
		}
		fprintf(fout, "\t\t\\node[anchor=south west] at (%f,%f) {\\csname %s%d\\endcsname};\n", x, y,
				Global::type_of_states == CS_63 ? "cs" : "bs", s);
		x += width + sep;
	}
	fprintf(fout, "\t\\end{tikzpicture}\n}");
	fclose(fout);
	delete[] acids;
}

