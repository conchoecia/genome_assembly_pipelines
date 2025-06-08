{
    cigar=$6; len=0;
    while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
        n=substr(cigar, RSTART, RLENGTH);
        val=substr(n, 1, length(n)-1) + 0;
        op=substr(n, length(n), 1);
        if (op ~ /[M=X]/) len += val;
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    if (len > 0) print $1, len
}