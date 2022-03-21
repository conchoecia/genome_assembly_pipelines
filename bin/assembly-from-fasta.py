
import os
import re
import sys
import pysam
import getopt

__authors__  = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Create an JBAT assembly file from a fasta file'


HELP_FLAGS  = set(('-h', '--help'))
MING_FLAGS  = set(('-g', '--min-gap-length'))
MAXG_FLAGS  = set(('-G', '--max-gap=length'))
MOVE_DEBRIS = set(('-D','--mv-Ns-to-debris'))
BROC_FLAGS  = set(('-c','--break-on-center'))

num = len

def _fold(linear_string, width=None):
    if width is None:
        return linear_string

    folded_string = []
    linear_length = len(linear_string)
    for offset in range(0, linear_length, width):
        if offset+width < linear_length:
            folded_string.append(linear_string[offset:offset+width])
        else:
            folded_string.append(linear_string[offset:])
            
    return '\n'.join(folded_string)


def format_fasta(name, sequence, width=None):
    return ">%s\n%s\n" % (name, _fold(sequence, width=width))


def write_record(record, contig_beg, contig_end, indices, debris=False, width=None, fasta=sys.stdout, assembly=sys.stderr):
    contig_name = "%s:::fragment_%d%s" % (record.name, indices[0], ':::debris' if debris else '')
    fasta.write(format_fasta(
        contig_name,
        record.sequence[contig_beg:contig_end],
        width=width
    ))
    assembly.write(">%s %d %d\n" % (contig_name, indices[1], contig_end-contig_beg))
        

def write_asm(scaffold_indices, assembly=sys.stderr):
    for scaffold in scaffold_indices:
        sep = ''
        for i in scaffold:
            assembly.write("%s%d" % (sep, i))
            sep = ' '
        assembly.write("\n")

    
def break_fasta(ifile, ofile, afile, min_gap_length=100, max_gap_length=100, mv_Ns_to_debris=False, break_on_center=False, line_length=None):
    nucleotides = re.compile('[acgturyswkmbdhvACGTURYSWKMBDHV]+')

    indices = [1, 1]
    contig_index = 0
    cprops_index = 1
    cprops_indices = []
    debris_indices = []
    for scaffold in ifile:
        scaffold_length = len(scaffold.sequence)
        contig_beg = 0
        prev_beg = 0
        prev_end = None

        indices[contig_index] = 1
        
        scaffold_indices = []
        cprops_indices.append(scaffold_indices)
        for contig in nucleotides.finditer(scaffold.sequence):
            curr_beg = contig.start()
            curr_end = contig.end()
            
            if prev_end is not None:
                gap_length = curr_beg - prev_end
                if min_gap_length <= gap_length <= max_gap_length:
                    if break_on_center:
                        curr_beg = prev_end + int(0.5 * gap_length)
                        write_record(
                            scaffold,
                            contig_beg,
                            curr_beg,
                            indices,
                            debris=False,
                            width=line_length,
                            fasta=ofile,
                            assembly=afile
                        )
                        scaffold_indices.append(indices[cprops_index])
                        indices[contig_index] += 1
                        indices[cprops_index] += 1

                        
                    else:
                        write_record(
                            scaffold,
                            contig_beg,
                            prev_end,
                            indices,
                            debris=False,
                            width=line_length,
                            fasta=ofile,
                            assembly=afile
                        )
                        scaffold_indices.append(indices[cprops_index])
                        indices[contig_index] += 1
                        indices[cprops_index] += 1
                        
                        write_record(
                            scaffold,
                            prev_end,
                            curr_beg,
                            indices,
                            debris=True,
                            width=line_length,
                            fasta=ofile,
                            assembly=afile
                        )
                        if mv_Ns_to_debris:
                            debris_indices.append(indices[cprops_index])
                        else:
                            scaffold_indices.append(indices[cprops_index])
                        
                        indices[contig_index] += 1
                        indices[cprops_index] += 1
                    
                    contig_beg = curr_beg
                    
            prev_beg = curr_beg
            prev_end = curr_end

        write_record(
            scaffold,
            contig_beg,
            scaffold_length,
            indices,
            debris=False,
            width=line_length,
            fasta=ofile,
            assembly=afile
        )
        scaffold_indices.append(indices[cprops_index])
        indices[contig_index] += 1
        indices[cprops_index] += 1
        
    cprops_indices.append(debris_indices)

    write_asm(cprops_indices, assembly=afile)


        
def usage(message=None, status=1):
    message = '' if message is None else 'ERROR: %s\n\n' % message
    sys.stderr.write("\n")
    sys.stderr.write("Program: %s (%s)\n" % (__program__, __purpose__))
    sys.stderr.write("Version: %s %s\n" % (__pkgname__, __version__))
    sys.stderr.write("Contact: %s\n" % __contact__)
    sys.stderr.write("\n")    
    sys.stderr.write("Description:\n\n")    
    sys.stderr.write("   Break a multi-fasta file of scaffolds into a fasta file of contigs,\n")
    sys.stderr.write("   also writes a JuiceBox assembly file describing the changes.\n")
    sys.stderr.write("\n")
    sys.stderr.write("Usage:   %s [options] <in.fasta> <out-prefix>\n" % (__program__))
    sys.stderr.write("\n")
    sys.stderr.write("Options: -c,--break-on-center       Break gaps on their center bases\n")
    sys.stderr.write("         -D,--mv-Ns-to-debris       Move gap N sequences to debris pile\n")
    sys.stderr.write("         -g,--min-gap-length <int>  Min gap size to break contigs [1]\n")
    sys.stderr.write("         -G,--max-gap-length <int>  Max gap size to break contigs [-1]\n")
    sys.stderr.write("         -h,--help                  Print this help message and exit\n")
    sys.stderr.write("\n")
    sys.stderr.write("Notes:  Set --max-gap-length to a value < 0 to break every gap larger than\n")
    sys.stderr.write("        --min-gap-length\n\n")
    sys.stderr.write("%s\n" % message)
    sys.exit(status)


def main(argv):
    try:
        options, arguments = getopt.getopt(argv, 'Dchg:G:',['mv-Ns-to-debris', 'help', 'min-gap-length', 'max-gap-length', 'break-on-center'])
    except getopt.GetoptError as message:
        usage(message)

    line_len = 100
    jat_compat = False
    min_gap_len = 1
    max_gap_len = -1
    mv_Ns_to_debris = False
    break_on_center = False
    for flag, value in options:
        try:
            if   flag in HELP_FLAGS: usage()
            elif flag in MING_FLAGS: min_gap_len = int(value)
            elif flag in MAXG_FLAGS: max_gap_len = int(value)
            elif flag in BROC_FLAGS: break_on_center = True
            elif flag in MOVE_DEBRIS: mv_Ns_to_debris = True
        except ValueError:
            usage("Invalid integer (via %r)" % flag)
        
    if num(arguments) < 2:
        usage('Too few arguments')
    if num(arguments) > 2:
        usage('Too many arguments')

    if max_gap_len < 0:
        max_gap_len = sys.maxsize
        
    if min_gap_len > max_gap_len:
        usage("{--max-gap-length} must be greater than or equal to {--min-gap-length}")
        
    ifile = pysam.FastxFile(arguments[0])
    ofile = open(arguments[1]+'.fasta', 'w')
    afile = open(arguments[1]+'.assembly', 'w')
    
    break_fasta(
        ifile,
        ofile,
        afile,
        min_gap_len,
        max_gap_len,
        mv_Ns_to_debris,
        break_on_center,
        line_len
    )

    ofile.close()
    afile.close()
    

if __name__ == '__main__':
    main(sys.argv[1:])
    
