from __future__ import print_function

import argparse
import collections
import gzip
import re
import sys


#from itertools import islice

def count_seqs(input_file, read_max_lines=None, top_hits=25, length=None):

    # This breaks read_max_lines but I don't care at the moment
    # that's for future Mario to worry about
    # suckerrr
    fparser = FastQParser(input_file)

    #with open(input_file, 'r') as f:
    #    for num, line in enumerate(f):
    #        if line.startswith('@'):
    #            first_headerline = num
    #            break

    #total_counts = 0
    #with open(input_file, 'r') as f:
    #    print("Reading input file \"{}\"...\n\n".format(input_file), file=sys.stderr)
    #    seq_counts = collections.Counter( (x[:length] for x in islice(f, first_headerline+1, read_max_lines, 4)) )
    seq_counts = collections.Counter([ read[1][:length] for read in fparser])
    total_counts = reduce( lambda x,y: x+y, seq_counts.values() )
        #for seq in islice(f, first_headerline, max_lines, 4):
        #    lengths_dict[ seq ] += 1

    len_longest_seq = len(max( seq_counts.keys(), key=len))


    str_formatter = "{1:<{0}}{2:>7}{3:>10}"
    num_formatter = "{1:<{0}}{2:>7}{3:>10.1%}"
    print("Total number of reads is {}".format(total_counts))
    print(str_formatter.format(len_longest_seq+5, "Sequence", "Counts", "Percent"))
    print(str_formatter.format(len_longest_seq+5, "-"*len_longest_seq, "-"*7, "-"*7))
    for seq, count in seq_counts.most_common(top_hits):
        print(num_formatter.format(len_longest_seq+5, seq.strip(), count, float(count)/total_counts))



class FastQParser:
    """Parser for fastq files, possibly compressed with gzip
       Iterates over one record at a time. A record consists
       of a list with 4 elements corresponding to 1) Header,
       2) Nucleotide sequence, 3) Optional header, 4) Qualities"""

    def __init__(self,file,filter=None):
        self.fname = file
        self.filter = filter
        fh = open(file,"rb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:
            self._fh = fh
        self._records_read = 0
        self._next = self.setup_next()

    def __iter__(self):
        return self

    def next(self):
        return self._next(self)

    def setup_next(self):
        """Return the function to return the next record
        """
        if self.filter is None or len(self.filter.keys()) == 0:
            def _next(self):
                self._records_read += 1
                return [self._fh.next().strip() for n in range(4)]
        else:
            def _next(self):
                while True:
                    record = [self._fh.next().strip() for n in range(4)]
                    header = parse_header(record[0])
                    skip = False
                    for k, v in self.filter.items():
                        if k in header and header[k] not in v:
                            skip = True
                            break
                    if not skip:
                        self._records_read += 1
                        return record
        return _next

    def name(self):
        return self.fname

    def rread(self):
        return self._records_read

    def seek(self,offset,whence=None):
        self._fh.seek(offset,whence)

    def close(self):
        self._fh.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser("Get the frequencies of the most common sequences appearing in a fastq file.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="The fastq file from which to read.")
    parser.add_argument("-m", "--read-max-lines", type=int, help="The maximum number of lines to read from the file. Defaults to entire file.")
    parser.add_argument("-t", "--top-hits", type=int, default=25, help="The number of hits to return. Defaults to 25.")
    parser.add_argument("-l", "--length", type=int, help="How many nts to count. Default entire read.")

    args = vars(parser.parse_args())
    count_seqs(**args)
