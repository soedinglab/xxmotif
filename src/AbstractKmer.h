#ifndef SRC_ABSTRACTKMER_H_
#define SRC_ABSTRACTKMER_H_

#include <stdint.h>

#include <functional>

#include <cassert>
#include <cstring>
#include <limits>
#include <string>
#include <sstream>

#include "utils.h"

class AbstractKmer {
    private:
        void operator=(const AbstractKmer&);

    public:
        struct id_type {
                bool isNumeric;
                union {
                        uint64_t num;
                        char* str;
                } val;
                id_type() {
                    isNumeric = true;
                }
                ~id_type() {
                    if (!isNumeric) {
                        delete[] val.str;
                    }
                }
                explicit id_type(const char* s) {
                    isNumeric = false;
                    val.str = new char[strlen(s) + 1];
                    strcpy(val.str, s);
                }
                id_type(const id_type &other) {
                    isNumeric = other.isNumeric;
                    if (isNumeric) {
                        val.num = other.val.num;
                    } else {
                        val.str = new char[strlen(other.val.str) + 1];
                        strcpy(val.str, other.val.str);
                    }
                }
                id_type& operator=(const id_type& other) {
                    if (!isNumeric) {
                        delete[] val.str;
                    }
                    isNumeric = other.isNumeric;
                    if (isNumeric) {
                        val.num = other.val.num;
                    } else {
                        val.str = new char[strlen(other.val.str) + 1];
                        strcpy(val.str, other.val.str);
                    }
                    return *this;
                }
                bool operator==(const id_type& other) const {
                    return (isNumeric && other.isNumeric && val.num
                            == other.val.num) || (!isNumeric
                            && !other.isNumeric && strcmp(val.str,
                            other.val.str) == 0);
                }
                bool operator<(const id_type& other) const {
                    if (isNumeric) {
                        if (other.isNumeric) {
                            return val.num < other.val.num;
                        } else {
                            return true;
                        }
                    } else {
                        if (other.isNumeric) {
                            return false;
                        } else {
                            return strcmp(val.str, other.val.str) < 0;
                        }
                    }
                }
                std::string toString() {
                    std::stringstream s;
                    if (isNumeric) {
                        s << val.num;
                    } else {
                        s << std::string(val.str);
                    }
                    return s.str();
                }
        };

    protected:
        id_type id;
        int _length;
        int _numMatches;

    public:
        virtual ~AbstractKmer() {
        }

        id_type getId() const {
            return id;
        }

        int length() const {
            return _length;
        }

        int numMatches() const {
            return _numMatches;
        }

        /**
         * Returns the nth state index, not counting wildcards.
         */
        virtual uint8_t charAt(const int n) const = 0;

        /**
         * Returns the number of wildcards after the given character.
         */
        virtual int gapsAfter(const int n) const = 0;

        /**
         *  Returns a list with the number of match positions in the motif
         */
        virtual motif_columns_type getMotifColumns() const = 0;

        /**
         * Changes the character at position offset to the given state.
         * The offset can be negative (grow to the left) or greater than
         * the length of the kmer (grow to the right).
         * The first character of the kmer has offset 0.
         */
        virtual void mutate(const int offset, const uint8_t state) = 0;

        virtual AbstractKmer* clone() const = 0;

        bool isMatchPosition(const int n) const;


        /**
         * String representation. States are represented by their index and
         * numBestCols profile cells.
         * I.e.: [01:F92|L02|00]XX[17:R86|K04|S01]XXXXX[31:S86|T02|A01]
         */
        virtual std::string toString(const int numBestCols = 0,
                const char* alphabet = NULL) const;

        /**
         * String representation, reducing states to their highest scoring
         * amino acid.
         * I.e.: FXXRXXXXXS
         */
        std::string bestAAString() const;
        std::string bestNucString() const;

        bool operator==(const AbstractKmer &other) const;

        static bool haveIdenticalGapPattern(const AbstractKmer &k1,
                const AbstractKmer &k2) {
            if (k1.length() != k2.length() || k1.numMatches()
                    != k2.numMatches()) {
                return false;
            } else {
                for (int i = 0; i < k1.numMatches(); ++i) {
                    if (k1.gapsAfter(i) != k2.gapsAfter(i)) {
                        return false;
                    }
                }
            }
            return true;
        }
};

class KmerPtrOrderCmp {
    public:
        bool operator()(const AbstractKmer *k1, const AbstractKmer *k2) {
            return k1->getId() < k2->getId();
        }
};

namespace std {
template<> struct hash<AbstractKmer::id_type> {
        size_t operator()(const AbstractKmer::id_type& id) const {
            if (id.isNumeric) {
                return std::hash<uint64_t>()(id.val.num);
            } else {
                /* We have to compute a hash value for the id string,
                 * not just for the pointer value.
                 * The following is taken from boost/functional/hash.cpp.
                 */
                const char *last = id.val.str + strlen(id.val.str) - 1;
                size_t seed = 0;
                for (const char *c = id.val.str; c<=last; ++c) {
                  seed ^= std::hash<char>()(*c) + 0x9e3779b9 + (seed<<6) + (seed>>2);
                }
                return seed;
            }
        }
};
}

std::ostream& operator<<(std::ostream& os, const AbstractKmer &k);
std::ostream& operator<<(std::ostream& os, const AbstractKmer::id_type &id);

#endif  // SRC_ABSTRACTKMER_H_
