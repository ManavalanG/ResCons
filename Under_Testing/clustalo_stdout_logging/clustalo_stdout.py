from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import subprocess, sys
import threading


cline = ClustalOmegaCommandline("clustalo", infile="seqs_more.fasta", iterations=3, verbose=True, force=True, outfile="out.fasta")

p = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=(sys.platform!="win32"))
while True:
    line = p.stdout.readline()
    if not line:
        break
    print line.strip('\n')

