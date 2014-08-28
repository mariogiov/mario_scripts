"""This isn't my finest work but I'm in a hurry. Anyway only God can judge me, so beat it.

Generate read counts using this bash script from the base dir:

for file in P1251_*/02_process-radtags/*fastq*; do num_lines=$(wc -l $file | cut -d " " -f 1) && num_reads=$(echo "$num_lines / 4" | bc) && echo -e "$num_reads\t$file"; done | tee -a 02_num_reads.tsv

"""
import collections
import re
import sys


old_sample_name = None
samples_dict = collections.defaultdict(list)
with open(sys.argv[1]) as fh:
    try:
        while fh:
            retained_list, discarded_list = fh.readline().split(), fh.readline().split()
            retained_quantity = int(retained_list[0])
            sample_name, lane_num, read_num = re.match(r'(P\d+_\d+)/\S*/(\d)_\S+_(\d)-untrimmed-retained.fastq', retained_list[1]).groups()
            discarded_quantity = int(discarded_list[0])

            samples_dict[sample_name].append({"lane_num": lane_num, "read_num": read_num,
                                              "retained_quantity": retained_quantity,
                                              "discarded_quantity": discarded_quantity})
    except IndexError:
        pass


for sample_name, sample_list in sorted(samples_dict.items()):
    print "====== Sample {} ======".format(sample_name)
    print "{:>3}\t{:>3}\t{:>10}\t{:>10}\t{:>10}\t{:>10}".format("read", "lane", "discarded", "retained", "total", "% discarded")
    total_sample_retained = 0
    total_sample_discarded = 0
    total_sample_reads = 0
    for fastq in sorted(sample_list, key=lambda x: x['read_num']):
        total_sample_retained += fastq["retained_quantity"]
        total_sample_discarded +=fastq["discarded_quantity"]
        total_quantity = fastq["discarded_quantity"] + fastq["retained_quantity"]
        total_sample_reads += total_quantity
        discarded_percentage = 100 * float(fastq["discarded_quantity"]) / total_quantity
        print "{:>3}\t{:>3}\t{:>10,}\t{:>10,}\t{:>10,}\t{:>10.2f}".format(fastq["read_num"], fastq["lane_num"],
                                                                      fastq["discarded_quantity"], fastq["retained_quantity"],
                                                                      total_quantity, discarded_percentage)
    total_sample_discarded_percentage = 100 * float(total_sample_discarded) / total_sample_reads
    print "\n{:^16}{:>10,}\t{:>10,}\t{:>10,}\t{:>10.2f}\n".format("total",
                                                                        total_sample_discarded,
                                                                        total_sample_retained,
                                                                        total_sample_reads,
                                                                        total_sample_discarded_percentage)
