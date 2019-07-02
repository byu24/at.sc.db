#!/global/projectb/scratch/bjcole/env_STARsolo/bin/python3
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse, sys
import itertools
import re

parser=argparse.ArgumentParser()
parser.add_argument('--in1', help='Barcode read file name')
parser.add_argument('--in2', help='Genomic read file name')
parser.add_argument('--run', help='Run name')
parser.add_argument('--out', help='Output file name')
args=parser.parse_args()


file_f = args.in1
file_r = args.in2
run = args.run
file_out = args.out

file_out_fh = open(file_out, "w")
count = 0

f_iter = FastqGeneralIterator(open(file_f,"rU"))
r_iter = FastqGeneralIterator(open(file_r,"rU"))

file_out_fh.write("@RG\tID:A\tSM:%s\n" % run)

for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter,r_iter):
    f_id_re = re.search('(\S+)', f_id)
    r_id_re = re.search('(\S+)', f_id)

    f_name = f_id_re.group(1)
    r_name = r_id_re.group(1)

    assert f_name == r_name
    cb = f_seq[0:12]
    um = f_seq[12:20]
    file_out_fh.write("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tXC:Z:%s\tRG:Z:A\tXM:Z:%s\n" % (r_name, r_seq, r_q, cb, um))
    count += 2

file_out_fh.close()
