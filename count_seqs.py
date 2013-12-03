import argparse
import collections
from itertools import islice
import re
import sys

def count_seqs(input_file, read_max_lines=None, top_hits=25):

    with open(input_file, 'r') as f:
        for num, line in enumerate(f):
            if line.startswith('@'):
                first_headerline = num
                break

    total_counts = 0
    with open(input_file, 'r') as f:
        sys.stderr.write("Reading input file \"{}\"...\n\n".format(input_file))
        seq_counts = collections.Counter( islice(f, first_headerline+1, read_max_lines, 4) )
        total_counts = reduce( lambda x,y: x+y, seq_counts.values() )
        #for seq in islice(f, first_headerline, max_lines, 4):
        #    lengths_dict[ seq ] += 1

    len_longest_seq = len( max( seq_counts.keys(), key=len) )

    str_formatter = "{1:<{0}}{2:>7}{3:>10}"
    num_formatter = "{1:<{0}}{2:>7}{3:>10.1%}"
    print( str_formatter.format(len_longest_seq+5, "Sequence", "Counts", "Percent"))
    print( str_formatter.format(len_longest_seq+5, "-"*len_longest_seq, "-"*7, "-"*7))
    for seq, count in seq_counts.most_common(top_hits):
        print num_formatter.format(len_longest_seq+5, seq.strip(), count, float(count)/total_counts)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Get the frequencies of the most common sequences appearing in a fastq file.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="The fastq file from which to read.")
    parser.add_argument("-m", "--read-max-lines", dest="read_max_lines", type=int, help="The maximum number of lines to read from the file. Defaults to entire file.")
    parser.add_argument("-t", "--top-hits", dest="top_hits", type=int, help="The number of hits to return. Defaults to 25.")

    args = vars(parser.parse_args())
    count_seqs(**args)
