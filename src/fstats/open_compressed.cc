#include "open_compressed.h"
#include "refcount_array.h"	/* refcount_array<> */
#include <cassert>	/* assert() */
#include <errno.h>	/* EINVAL, EISDIR, ENFILE, ENOENT, errno */
#include <fcntl.h>	/* O_RDONLY, open() */
#include <list>		/* list<> */
#include <map>		/* map<> */
#include <stdio.h>	/* fprintf(), stderr */
#include <stdlib.h>	/* NULL */
#include <string.h>	/* strchr(), strerror() */
#include <string>	/* string */
#include <sys/stat.h>	/* S_IFDIR, stat(), struct stat */
#include <sys/types.h>	/* pid_t, size_t, ssize_t */
#include <sys/wait.h>	/* WNOHANG, waitpid() */
#include <unistd.h>	/* _SC_OPEN_MAX, _exit(), close(), dup2(), execlp(), fork(), pipe(), read(), sysconf() */

#define BUFSIZE 32768

// all function calls have to be inlined to prevent conflicts with
// write_fork, which also uses a LocalData class

class LocalData {
    private:
	/* map of open files to gzip process id's */
	std::map<int, pid_t> open_processes;
	/* list of closed processes that need to be waited on */
	std::list<pid_t> closed_processes;
	void finish_nohang(void) {
		std::list<pid_t>::iterator __a = closed_processes.begin();
		const std::list<pid_t>::iterator __end_a = closed_processes.end();
		while (__a != __end_a) {
			if (waitpid(*__a, NULL, WNOHANG) > 0) {
				__a = closed_processes.erase(__a);
			} else {
				++__a;
			}
		}
	}
    public:
	const long open_max;
	/* internal file buffers - has to be a raw array, sadly */
	refcount_array<char> *buffers;
	ssize_t *buffer_start;
	ssize_t *buffer_length;
	LocalData(void) : open_max(sysconf(_SC_OPEN_MAX)) {
		assert(open_max != 0);
		buffers = new refcount_array<char>[open_max];
		assert(buffers != NULL);
		buffer_start = new ssize_t[open_max];
		assert(buffer_start != NULL);
		buffer_length = new ssize_t[open_max];
		assert(buffer_length != NULL);
	}
	void add_open(const int __i, const pid_t __j) {
		open_processes[__i] = __j;
	}
	void close_process(const int __i) {
		assert(-1 < __i && __i < open_max);
		close(__i);
		buffers[__i].resize(0);
		std::map<int, pid_t>::iterator __a = open_processes.find(__i);
		if (__a != open_processes.end()) {
			closed_processes.push_back(__a->second);
			open_processes.erase(__a);
		}
		finish_nohang();
	}
	~LocalData(void) {
		delete[] buffers;
		delete[] buffer_start;
		delete[] buffer_length;
		std::map<int, pid_t>::const_iterator __a = open_processes.begin();
		const std::map<int, pid_t>::const_iterator __end_a = open_processes.end();
		for (; __a != __end_a; ++__a) {
			close(__a->first);
		}
		std::list<pid_t>::const_iterator __b = closed_processes.begin();
		const std::list<pid_t>::const_iterator __end_b = closed_processes.end();
		for (; __b != __end_b; ++__b) {
			waitpid(*__b, NULL, 0);
		}
		for (__a = open_processes.begin(); __a != __end_a; ++__a) {
			waitpid(__a->second, NULL, 0);
		}
	}
};

static LocalData local;

// simply return the suffix of the file name, if it matches one of the list

void get_suffix(const std::string &filename, std::string &suffix) {
	const char *suffix_list[3] = { ".gz", ".bz2", ".Z" };
	const size_t suffix_size[3] = { 3, 4, 2 };
	suffix.clear();
	for (int i = 0; i != 3; ++i) {
		const size_t &j = suffix_size[i];
		if (filename.size() > j && filename.compare(filename.size() - j, j, suffix_list[i]) == 0) {
			suffix = suffix_list[i];
			break;
		}
	}
}

/*
 * returns empty, .Z, .gz, or .bz2; checks to see if the filename ends in any
 * given suffix, if not checks to see if a file with the given suffix exists
 */

int find_suffix(std::string &filename, std::string &suffix) {
	const char *suffix_list[3] = { ".gz", ".bz2", ".Z" };
	const size_t suffix_size[3] = { 3, 4, 2 };
	suffix.clear();
	for (int i = 0; i != 3; ++i) {
		const size_t &j = suffix_size[i];
		if (filename.size() > j && filename.compare(filename.size() - j, j, suffix_list[i]) == 0) {
			suffix = suffix_list[i];
			break;
		}
	}
	struct stat buf;
	if (stat(filename.c_str(), &buf) == 0) {
		/* only open regular files */
		if ((buf.st_mode & S_IFDIR) != 0) {
			errno = EISDIR;
			return -1;
		} else {
			return 0;
		}
	} else if (errno != ENOENT) {
		fprintf(stderr, "Error: stat: %s: %s\n", filename.c_str(), strerror(errno));
		return -1;
	} else if (!suffix.empty()) {
		return -1;
	}
	for (int i = 0; i != 3; ++i) {
		const std::string s(filename + suffix_list[i]);
		if (stat(s.c_str(), &buf) == -1) {
			if (errno != ENOENT) {
				fprintf(stderr, "Error: stat: %s: %s\n", filename.c_str(), strerror(errno));
				return -1;
			}
		} else if ((buf.st_mode & S_IFDIR) == 0) {
			suffix = suffix_list[i];
			filename += suffix;
			return 0;
		}
	}
	errno = ENOENT;
	return -1;
}

/*
 * open a file, with gzip if the filename ends in .gz, .bz2, or .Z; if file is
 * not found, .gz, .bz2, and .Z are added to end, to see if file is compressed
 */

int open_compressed(const std::string &filename) {
	int fd = -1;
	std::string s(filename);
	std::string suffix;
	/* see if file exists */
	if (!s.empty() && s.compare("-") != 0 && find_suffix(s, suffix) == -1) {
		return -1;
	}
	if (!suffix.empty()) {
		pid_t pid;
		int pipefd[2];
		if (pipe(pipefd) == -1) {
			fprintf(stderr, "Error: pipe: %s\n", strerror(errno));
			return -1;
		} else if (pipefd[0] >= local.open_max) { // too many open files
			close(pipefd[0]);
			close(pipefd[1]);
			errno = ENFILE;
			fprintf(stderr, "Error: open: %s\n", strerror(errno));
			return -1;
		} else if ((pid = fork()) == -1) {
			fprintf(stderr, "Error: fork: %s\n", strerror(errno));
			close(pipefd[0]);
			close(pipefd[1]);
			return -1;
		} else if (pid == 0) {	/* child */
			close(0);
			close(pipefd[0]);
			if (dup2(pipefd[1], 1) == -1) {
				fprintf(stderr, "Error: dup2: %s\n", strerror(errno));
				_exit(1);
			}
			close(pipefd[1]);
			for (int i = 3; i < local.open_max; ++i) {
				close(i);
			}
			if (suffix == ".bz2") {
				if (execlp("bzip2", "bzip2", "-d", "-c", s.c_str(), NULL) == -1) {
					fprintf(stderr, "Error: execlp bzip2 -d -c: %s\n", strerror(errno));
				}
			} else {
				if (execlp("gzip", "gzip", "-d", "-c", s.c_str(), NULL) == -1) {
					fprintf(stderr, "Error: execlp gzip -d -c: %s\n", strerror(errno));
				}
			}
			_exit(1);
		} else {		/* parent */
			close(pipefd[1]);
			fd = pipefd[0];
			local.add_open(fd, pid);
		}
	} else if (s.empty() || s.compare("-") == 0) {
		fd = 0;
	} else if ((fd = open(s.c_str(), O_RDONLY)) == -1) {
		fprintf(stderr, "Error: open: %s\n", strerror(errno));
		return -1;
	} else if (fd >= local.open_max) {	// too many open files
		close(fd);
		errno = ENFILE;
		fprintf(stderr, "Error: open: %s\n", strerror(errno));
		return -1;
	}
	// make sure to leave extra spot for trailing null
	local.buffers[fd].resize(BUFSIZE + 1);
	local.buffers[fd][0] = '\0';
	local.buffer_start[fd] = 0;
	local.buffer_length[fd] = 0;
	return fd;
}

/* close the file and wait on the gzip process, if any */

void close_compressed(const int fd) {
	local.close_process(fd);
}

/*
 * read input from a file descriptor - returns data up to end of line or
 * end of file (and strips the end of line character); returns -1 on error,
 * otherwise number of characters read (minus end of line, if any)
 */

ssize_t pfgets(const int fd, std::string &line) {
	if (fd < 0 || local.open_max <= fd) {
		fprintf(stderr, "Error: pfgets: fd out of range: %d\n", fd);
		return -1;
	}
	if (local.buffers[fd].empty()) {
		fprintf(stderr, "Error: pfgets: buffer unallocated\n");
		return -1;
	}
	line.clear();
	char * const buf = local.buffers[fd].array();
	ssize_t &i = local.buffer_start[fd];
	ssize_t &j = local.buffer_length[fd];
	for (;;) {
		char * const end_of_line = strchr(buf + i, '\n');
		if (end_of_line != NULL) {
			*end_of_line = '\0';
			line += buf + i;
			i = end_of_line - buf + 1;
			return (ssize_t)line.size();
		}
		line += buf + i;
		i = 0;
		j = read(fd, buf, BUFSIZE);
		if (j <= 0) {
			if (j == -1) {
				fprintf(stderr, "Error: read(%d): %s\n", fd, strerror(errno));
			}
			buf[j = 0] = '\0';
			/* only return -1 if we don't return anything else */
			return line.empty() ? -1 : (ssize_t)line.size();
		}
		buf[j] = '\0';
	}
}

/* read up to size bytes from fd and put them into ptr */

ssize_t pfread(const int fd, void *ptr, const size_t size) {
	if (fd < 0 || local.open_max <= fd) {
		fprintf(stderr, "Error: pfgets: fd out of range: %d\n", fd);
		return -1;
	}
	if (local.buffers[fd].empty()) {
		fprintf(stderr, "Error: pfgets: buffer unallocated\n");
		return -1;
	}
	char *s = (char *)ptr;
	size_t k = size;
	char * const buf = local.buffers[fd].array();
	ssize_t &i = local.buffer_start[fd];
	ssize_t &j = local.buffer_length[fd];
	for (;;) {
		if ((size_t)(j - i) >= k) {
			memcpy(s, buf + i, k);
			i += k;
			return size;
		}
		const size_t n = j - i;
		memcpy(s, buf + i, n);
		s += n;
		k -= n;
		i = 0;
		j = read(fd, buf, BUFSIZE);
		if (j <= 0) {
			if (j == -1) {
				fprintf(stderr, "Error: read(%d): %s\n", fd, strerror(errno));
			}
			buf[j = 0] = '\0';
			/* only return -1 if we don't return anything else */
			return k == size ? -1 : (ssize_t)(size - k);
		}
		buf[j] = '\0';
	}
}

// like pfread(), but doesn't use buffer; obviously, don't use along with
// pfread(); currently intended for testing purposes only

ssize_t pfread2(const int fd, void *ptr, const size_t size) {
	assert(-1 < fd && fd < local.open_max);
	char *buf = (char *)ptr;
	ssize_t i = size;
	while (i != 0) {
		const ssize_t j = read(fd, buf, i);
		if (j == -1) {
			fprintf(stderr, "Error: read (%d %p): %s\n", fd, ptr, strerror(errno));
			return -1;
		}
		i -= j;
		buf += j;
	}
	return size;
}
