#ifndef _OPEN_COMPRESSED_H
#define _OPEN_COMPRESSED_H

#include <string>	/* string */
#include <sys/types.h>	/* ssize_t */

extern void get_suffix(const std::string &, std::string &);
extern int find_suffix(std::string &, std::string &);
extern int open_compressed(const std::string &);
extern void close_compressed(const int);
extern ssize_t pfgets(const int, std::string &);
extern ssize_t pfread(const int, void *, const size_t);
extern ssize_t pfread2(const int, void *, const size_t);

#endif /* !_OPEN_COMPRESSED_H */
