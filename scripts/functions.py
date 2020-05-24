"""
these are functions that are common to several pipelines
"""
from subprocess import Popen, PIPE
import subprocess

def percent_of_data_mapping(bam, read_id_length, output_txt):
    """
    This function looks in a bam file calculates the percent
      of data in the txt file that map to that bam file.
    """
    transcript_to_len = {}
    added = set()
    tot_mapped = 0
    num_reads = 0
    tot_size = 0
    print("reading in the read id to read length file")
    with open(read_id_length, "r") as f:
        for line in f:
            if line.strip():
                splitd = line.split()
                transcript_to_len[splitd[0]] = int(splitd[1])
                tot_size += int(splitd[1])
                num_reads += 1
    print("done reading in the read ids to lengths")
    command = "samtools view -F 4 {} | cut -f1".format(bam)
    process = Popen(command, stdout=PIPE, shell=True)
    blank_counter = 0
    while True:
        line = process.stdout.readline().rstrip().decode("utf-8").strip()
        if (len(added) % 1000) == 0:
           sys.stdout.write("  {0:.3f}% Done. {1:.3f}% of data mapped.\r".format(
                            100 * (len(added)/num_reads),
                            100 * (tot_mapped/tot_size)))
           sys.stdout.flush()
        if not line.strip():
            blank_counter += 1
            if blank_counter > 99999:
                break
        else:
            if line not in added:
                added.add(line)
                tot_mapped += transcript_to_len[line]

    with open(output_txt, "w") as f:
        print("Number of bases mapped: {}".format(tot_mapped), file=f)
        print("Number of bases total: {}".format(tot_size), file=f)
        print("{0:.4}%".format(100* (tot_mapped/tot_size)), file=f)
        print("", file=f)
        print("Number of reads mapped: {}".format(len(added)), file=f)
        print("Number of reads total: {}".format(num_reads), file=f)
        print("{0:.4}%".format(100* (len(added)/num_reads)), file=f)

def run_fasta_stats(fpath):
    """
    This runs fasta_stats on an assembly and returns a dictionary of the values.
    """
    this_data = {}
    # This version is if you want N95_scaflen. Not a good idea for bad assemblies.
    #new_fields = ["num_scaffolds", "num_tigs", "tot_size_scaffolds",
    #              "tot_size_tigs", "scaffold_N50", "scaffold_L50",
    #              "contig_N50", "contig_L50", "perGap", "N95_scaflen", "numGaps"]
    new_fields = ["num_scaffolds",
                  "num_tigs",
                  "tot_size_scaffolds",
                  "tot_size_tigs",
                  "scaffold_N50",
                  "scaffold_L50",
                  "scaffold_N90",
                  "scaffold_L90",
                  "scaffold_N95",
                  "scaffold_L95",
                  "contig_N50",
                  "contig_L50",
                  "contig_N90",
                  "contig_L90",
                  "contig_N95",
                  "contig_L95",
                  "perGap",
                  "numGaps"]
    tcmd = "{} {}".format(fs_path, fpath).split(" ")
    results = subprocess.run(tcmd, stdout=subprocess.PIPE).stdout.decode('utf-8').split()
    assert len(new_fields)==len(results)
    for i in range(len(results)):
        if new_fields[i] == "perGap":
            this_data[new_fields[i]] = float(results[i].strip('%'))
        else:
            this_data[new_fields[i]] = int(results[i])
    return(this_data)

def parse_flagstat(tfile):
    """
    just gets the percent mapping from flagstat
    """
    pmap = -1
    with open(tfile, "r") as f:
        counter = 0
        for line in f:
            if line.strip():
                counter += 1
                if counter == 5:
                    assert "mapped" in line
                    pmap = float(line.split("(")[1].split("%")[0])
                    break
    return pmap

def parse_permap(tfile):
    """
    just gets the percent of bases in the reads that mapped,
      and the number of reads that mapped
    """
    bmap = -1
    pmap = -1
    with open(tfile, "r") as f:
        counter = 0
        for line in f:
            line = line.strip()
            if line:
                counter += 1
                if counter == 3:
                    assert "%" in line
                    bmap = float(line.replace("%",""))
                elif counter == 6:
                    assert "%" in line
                    pmap = float(line.replace("%",""))

    return [bmap, pmap]
