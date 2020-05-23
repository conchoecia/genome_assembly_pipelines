#include "itoa.h"	/* itoa() */
#include "pretty_print.h"
#include <string>	/* string */
#include <sys/types.h>	/* size_t */

/*
 * converts a number into a string (like itoa()), but it sticks in
 * commas every three numbers (making it easier to read big numbers)
 */

std::string pretty_print(unsigned long x, bool negative) {
	std::string s = itoa(x, negative);
	int commas = (s.length() - 1 - (negative ? 1 : 0)) / 3;
	size_t a = s.length() - 1;
	size_t b = a + commas;
	s.resize(b + 1);
	int i;
	for (i = 0; a < b; --a, --b, ++i) {
		s[b] = s[a];
		if (i == 2) {
			i = -1;
			s[--b] = ',';
		}
	}
	return s;
}
