#!/usr/bin/python3

import argparse

class Sequence:
    def __init__(self, name, line):
        self.name = name
        self.line = line
        self.seq = ''
        self.size = 0
        
    def set_sequence(self, seq, store=False, fast=False):
        if fast:
            self.seq += seq
            return
        if self.size < 1:
            self.seq += self.seq[:10]
        elif store:
            self.seq += seq
        self.size += len(seq)
        
    def get_name(self):
        return self.name
    
    def get_line(self):
        return self.line
    
    def get_seq(self):
        return self.seq
    

class Fasta:
    def __init__(self, arquivo):
        self.arquivo = arquivo
        self.sequences = []
        
    def get_seqs(self):
        return self.sequences
        
    def load_some(self, seqs):
        seq, cont, loaded, lazy = None, 0, [], len(seqs) > 0
        with open(self.arquivo) as fin:
            for ln in fin:
                cont += 1
                if ln.startswith('>'):
                    if lazy and all([x in loaded for x in seqs]):
                        return self.sequences
                    s = ln[1:].split(' ')[0]
                    seq = Sequence(s, cont)
                    self.sequences.append(seq)
                    loaded.append(s)
                elif seq:
                    seq.set_sequence(ln.strip(), lazy)
        return self.sequences
    
    def load(self):
        return self.load_some([])
        
def MAIN(file, seqs):
    fasta = Fasta(file)
    fasta.load_some(seqs)

    for s in fasta.sequences:
        print(s.name, s.line, s.size)

    fasta.load_some(fasta.sequences[0].name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Script to list or extract fasta records.
    output mode list:
        name,line,size
    output mode extract:
        name,line,size,sequence
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('fasta')
    parser.add_argument('-s', '--seqs', nargs='+', default=[])
    args = parser.parse_args()
    MAIN(args.fasta, args.seqs)