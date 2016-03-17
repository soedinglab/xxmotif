// Copyright 2009, Andreas Biegert

#ifndef SRC_LOG_LEVEL_H_
#define SRC_LOG_LEVEL_H_

#include <cstdio>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <string>

/* usage: -DLOG_MAX_LEVEL=2
				  0		  1		  2		3	   4	   5	   6	   7 */
enum LogLevel { ERROR, WARNING, INFO, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4 };

class Log {
 public:
  Log() {}
  virtual ~Log();
  std::ostringstream& get(LogLevel level = INFO);

  static LogLevel& reporting_level();
  static std::string to_string(LogLevel log_level);
  static LogLevel from_string(const std::string& log_level);
  static LogLevel from_int(int log_level);
  static FILE*& stream() {
    static FILE* p_stream = stderr;
    return p_stream;
  }

 protected:
  std::ostringstream os;
  LogLevel level;

 private:
  Log(const Log&);
  Log& operator =(const Log&);
};

#ifndef LOG_MAX_LEVEL
#define LOG_MAX_LEVEL DEBUG4
#endif

#define LOG(level)                                              \
  if (level > LOG_MAX_LEVEL || !kDebug) ;                       \
  else if (level > Log::reporting_level() || !Log::stream()) ;  \
  else Log().get(level)


#endif  // SRC_LOG_LEVEL_H_
