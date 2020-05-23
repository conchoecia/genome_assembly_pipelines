#include "itoa.h"
#include <math.h>	/* floor(), log10() */
#include <string>	/* string */

/* converts an unsigned long (with option sign) to a string */

std::string itoa(unsigned long x, bool negative) {
	int length = negative ? 2 : 1;
	if (x > 9) {
		length += (int)floor(log10(x));
	}
	std::string value;
	value.resize(length);
	std::string::reverse_iterator a = value.rbegin();
	std::string::reverse_iterator end = value.rend();
	if (negative) {
		value[0] = '-';
		--end;
	}
	for (; a != end; ++a, x /= 10) {
		*a = '0' + (x % 10);
	}
	return value;
}
