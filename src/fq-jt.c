#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <zlib.h>

#define MAX_ID_LEN 255
#define MAX_SEQ_LEN 1023
#define DEFAULT_CHUNK_SIZE 1024
#define OFFSET DEFAULT 0
/* Type to hold the forward and reverse read
   of a sequence pair with quality scores */
typedef struct sqp {
  char id[MAX_ID_LEN+1];
  char seq[MAX_SEQ_LEN+1];
  char qual[MAX_SEQ_LEN+1];
  char seq_trunc[MAX_SEQ_LEN+1];
  char qual_trunc[MAX_SEQ_LEN+1];
} Sqp;
typedef struct sqp* SQP;

/* Type to hold an array of sqp, i.e, a database
   of sequences pairs and their quality scores */
typedef struct sqpdb {
  SQP* sqps; // pointer to array of seq & qual pairs
  unsigned int seq_len; // the length of the sequences
  size_t num_reads; // The current length of the database, 0-indexed
  size_t size; // The currently allocated size of the database
} Sqpdb;
typedef struct sqpdb* SQPDB;

int is_gz( const char* fq_fn );
void write_trunc_seq( const SQP curr_seq );
int trunc_sqp( SQP seq_p, const char* junc_seq, const int offset );
SQPDB init_SQPDB( size_t size );
int read_fastq( FILE* fastq, SQP curr_seq );
int gz_read_fastq( gzFile fastq, SQP curr_seq );
FILE * fileOpen(const char *name, char access_mode[]);
void help( void ) {
  printf( "fq-jt -f <fastq file> -t <truncation sequence> -l <offset length\n" );
  printf( "Takes an input fastq file and removes any sequence after the\n" );
  printf( "truncation sequence, if it's present. If an offset length is\n" );
  printf( "specified, this much sequence is kept, past the start of the\n" );
  printf( "position of the truncation sequence.\n" );
  exit( 0 );
}

/*** OVERVIEW ***/
/* Read in CHUNK_SIZE fastq sequences into memory
   Set up CHUNK_SIZE number of threads, handing each a sequence to truncate
   Call pthread_join on head thread
   Write out the CHUNK_SIZE truncated fastq sequence
*/

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  char fq_fn[MAX_ID_LEN+1];
  char junc_seq[MAX_ID_LEN+1];
  int ich;
  int offset;
  FILE* fq_p;
  gzFile fq_gz_p;
  SQP curr_seq;

  if ( argc == 1 ) {
    help();
  }
  while( (ich=getopt( argc, argv, "f:t:l:h" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_fn, optarg );
      break;
    case 't' :
      strcpy( junc_seq, optarg );
      break;
    case 'l' :
      offset = atoi( optarg );
      break;
    case 'h' :
      help();
    }
  }

  /* Initialize the SQP curr_seq */
  curr_seq = (SQP)malloc(sizeof( Sqp ));

  /* Open the fastq file */
  if ( is_gz( fq_fn ) ) {
    fq_gz_p = gzopen( fq_fn, "r" );
    if ( fq_gz_p == NULL ) {
      help();
    }
    while( gz_read_fastq( fq_gz_p, curr_seq ) ) {
      trunc_sqp( curr_seq, junc_seq, offset );
      write_trunc_seq( curr_seq );
    }
    gzclose( fq_gz_p );
  }

  else {
    fq_p = fileOpen( fq_fn, "r" );
    if ( fq_p == NULL ) {
      help();
    }
    while( read_fastq( fq_p, curr_seq ) ) {
      trunc_sqp( curr_seq, junc_seq, offset );
      write_trunc_seq( curr_seq );
    }
    fclose( fq_p );
  }

  free( curr_seq );

  exit( 0 );
}

int is_gz( const char* fq_fn ) {
  size_t fn_len;
  fn_len = strlen( fq_fn );
  if ( (fq_fn[fn_len-3] == '.') &&
       (fq_fn[fn_len-2] == 'g') &&
       (fq_fn[fn_len-1] == 'z') ) {
    return 1;
  }
  return 0;
}

/* trunc_sqp
   ARGS: SQP seq_p - populated SQP with a sequence to be on
         the junction sequence if it occurs
     char* junc_seq - the junction sequence to look for
   RETURNS: 0 if no problems
   Populates the seq_trunc field of the input SQP
*/
int trunc_sqp( SQP seq_p, const char* junc_seq, const int offset ) {
  char* found;
  size_t truncated_len;
  found = strstr( seq_p->seq, junc_seq );
  if ( found == NULL ) {
    strcpy( seq_p->seq_trunc, seq_p->seq );
    strcpy( seq_p->qual_trunc, seq_p->qual );
  }
  else {
    truncated_len = found - &seq_p->seq[0];
    strncpy( seq_p->seq_trunc, seq_p->seq, truncated_len + offset );
    strncpy( seq_p->qual_trunc, seq_p->qual, truncated_len + offset);
    seq_p->seq_trunc[truncated_len + offset] = '\0';
    seq_p->qual_trunc[truncated_len + offset] = '\0';
  }
  return 0;
}

/* init_SQPDB
   Initialize and return a SQPDB
   Args: (1) size_t size - how big the array of sqps should be
*/
SQPDB init_SQPDB( size_t size ) {
  size_t i;
  SQPDB sqpdb;
  SQP first_seq;

  /* Try to allocate memory for SQPDB */
  sqpdb = (SQPDB)malloc(sizeof(Sqpdb));
  if ( sqpdb == NULL ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  }

  /* Try to allocate memory for sqps array */
  first_seq = (SQP)malloc( size * sizeof(Sqp) );
  sqpdb->sqps = (SQP*)malloc( size * sizeof(SQP) );
  
  if ( (first_seq   == NULL) ||
       (sqpdb->sqps == NULL) ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  } 
    
  for( i = 0; i < size; i++ ) {
    sqpdb->sqps[i] = &first_seq[i];
  }

  sqpdb->seq_len   = MAX_SEQ_LEN;
  sqpdb->num_reads = 0;
  sqpdb->size      = size;

  return sqpdb;
}

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( FILE* fastq, SQP curr_seq ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    curr_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      curr_seq->id[i] = '\0';
    }
  }
  curr_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      curr_seq->seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  curr_seq->seq[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", curr_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      curr_seq->qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  curr_seq->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return 0;
  }
  return 1;
}

/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int gz_read_fastq( gzFile fastq, SQP curr_seq ) {
  char c;
  size_t i;
  c = gzgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=gzgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    curr_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      curr_seq->id[i] = '\0';
    }
  }
  curr_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = gzgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = gzgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      curr_seq->seq[i++] = c;
    }
    c = gzgetc( fastq );
  }
  curr_seq->seq[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = gzgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", curr_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = gzgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = gzgetc( fastq );
  }

  /* Now, get the quality score line */
  c = gzgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      curr_seq->qual[i++] = c;
    }
    c = gzgetc( fastq );
  }
  curr_seq->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return 0;
  }
  return 1;
}


/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}

void write_trunc_seq( const SQP curr_seq ) {
  printf( "@%s\n", curr_seq->id );
  printf( "%s\n", curr_seq->seq_trunc );
  printf( "+\n" );
  printf( "%s\n", curr_seq->qual_trunc );
}
