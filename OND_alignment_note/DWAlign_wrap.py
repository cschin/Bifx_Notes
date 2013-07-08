from pbcore.io import FastaReader
from ctypes import *



class Alignment(Structure):
    """
    typedef struct {    
        seq_size aln_str_size ;
        seq_size dist ;
        char * t_aln_str;
        char * q_aln_str;
    } alignment;
    """
    _fields_ = [ ("aln_str_size", c_long),
                 ("dist", c_long),
                 ("t_aln_str", c_char_p),
                 ("q_aln_str", c_char_p)]


DWA = CDLL("./DWAlign.so")
DWA.align.argtypes = [ POINTER(c_char), c_long, POINTER(c_char), c_long, c_long, c_long ] 
DWA.align.restype = POINTER(Alignment)


seq = ""
f = FastaReader("test.fa")
for r in f:
    seq = r.sequence.upper() * 5

import random
random.seed(42)
sim_seq = []
ins_rate = 0.01
del_rate = 0.01
mis_rate = 0.01


for c in seq:
    while 1:
        if random.uniform(0, 1) < ins_rate:
            sim_seq.append( random.choice(["A","C","G","T"]))
        else:
            break
    if random.uniform(0, 1) < del_rate:
        continue
    if random.uniform(0, 1) < mis_rate:
        sim_seq.append( random.choice(["A","C","G","T"]))
        sim_seq.append(c)
        continue
    sim_seq.append(c)
sim_seq = "".join(sim_seq)

    
sim_seq2 = []
for c in seq:
    while 1:
        if random.uniform(0, 1) < ins_rate:
            sim_seq2.append( random.choice(["A","C","G","T"]))
        else:
            break
    if random.uniform(0, 1) < del_rate:
        continue
    if random.uniform(0, 1) < mis_rate:
        sim_seq2.append( random.choice(["A","C","G","T"]))
        sim_seq2.append(c)
        continue
    sim_seq2.append(c)
sim_seq2 = "".join(sim_seq2)


MAX_D = len(sim_seq) + len(sim_seq2)

for i in range(100):
    alignment = DWA.align(sim_seq, len(sim_seq),
                          sim_seq2, len(sim_seq2),
                          2000, 60)
    if i != 99:
        DWA.free_alignment(alignment)

print alignment[0].aln_str_size
print alignment[0].dist
print alignment[0].q_aln_str
print alignment[0].t_aln_str
