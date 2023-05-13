#!/usr/bin/python3

import argparse
from Fasta import Fasta


class GxF:

    def __init__(self, file):
        self.file = file
        error = False

        try:
            gtf, gff = False, False
            smp = [x for x in open(file).readlines(1000)
                   if x.startswith('#')][0]
            gtf = file.endswith('.gtf') and smp.index('#gtf') >= 0
            gff = (file.endswith('.gff') or file.endswith(
                '.gff3')) and smp.index('#gff') >= 0
        except:
            error = True
        if error or (not gtf and not gff):
            raise Exception('Impossivel inferir tipo GFF/GTF: ' + file)

        self.gff = gff
        self.gtf = not gff

    def convert(self, genes=[]):
        mrna2gene, genes, found = {}, set(genes), set()
        T = 'mRNA' if self.gff else 'transcript'

        def get(e, f):
            if f in e:
                return e.split(f)[1].split(';')[0].replace('"', '').strip()

        with open(self.file) as fin:
            for line in fin:
                if T in line:
                    fs = line.split(chr(9))
                    if len(fs) == 9 and fs[2] == T:
                        gene = get(
                            fs[-1], 'Parent=' if self.gff else 'gene_id ')
                        transcript = get(
                            fs[-1], 'ID=' if self.gff else 'transcript_id ')
                        if gene and transcript and gene in genes:
                            mrna2gene[transcript] = gene

        ncbi = all([x.startswith('rna-') for x in mrna2gene])

        with open(self.file) as fin, open(self.file+('.gtf' if self.gff else '.gff'), 'w') as fout:
            for line in fin:
                if len(line) < 8 or line[0] == '#':
                    continue
                f1, f2, f3, f4, f5, f6, f7, f8, f9 = line.split(chr(9))
                gene, transcript = None, None

                if f3 == 'gene':
                    gene = get(f9, 'ID=' if self.gff else 'gene_id ')
                elif f3 == T:
                    gene = get(f9, 'Parent=' if self.gff else 'gene_id ')
                    transcript = get(
                        f9, 'ID=' if self.gff else 'transcript_id ')
                elif f3 == 'exon' or f3 == 'CDS':
                    transcript = get(
                        f9, 'Parent=' if self.gff else 'transcript_id ')
                    gene = mrna2gene[transcript] if transcript in mrna2gene else None
                else:
                   continue

                if gene:
                    if gene in genes:
                        found.add(gene)
                    else:
                        continue
                    gene = gene[5:] if ncbi else gene
                else:
                    gene = ''

                if transcript:
                    if not transcript in mrna2gene:
                        continue
                    transcript = transcript[4:] if ncbi else transcript
                else:
                    transcript = ''

                fx = ''
                if self.gtf:
                    if f3 == 'gene':
                        fx = f'ID={gene};'
                    elif f3 == 'transcript':
                        f3 = 'mRNA'
                        fx = f'ID={transcript};Parent={gene};'
                    elif f3 == 'exon' or f3 == 'CDS':
                        fx = f'Parent={transcript};'

                fout.write(chr(9).join([
                    f1, f2, f3, f4, f5, f6, f7, f8,
                    f'gene_id "{gene}"; transcript_id "{transcript}";' if self.gff else fx
                ]) + '\n')
        return mrna2gene, found


def MAIN(file, w=None, x=False):
    gns = set()
    if w:
        with open(w) as fin:
            for l in fin:
                l = l.strip()
                if len(l) > 1:
                    gns.add(l)
            print(f'Filter for {len(gns)} genes ... ', end='')
    gxf = GxF(file)
    mrna2gene, found = gxf.convert(gns)
    if len(found) != len(gns):
        print(gns.difference(found), ' NOT FOUND.')
    if x:
        with open(x, 'w') as fo:
            for m, g in mrna2gene.items():
                fo.write(f'{g},{m}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Script clean, format and filter GFF files.
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('gxf', help="GTF or GFF file")
    # parser.add_argument('-f', '--fasta', help='Fasta with genomic seqs')
    # parser.add_argument('-b', '--black', help="black list to remove genes")
    parser.add_argument(
        '-w', '--white', help="white list to selecet this genes")
    # parser.add_argument('-g', '--gene', help='file for write genes seqs')
    # parser.add_argument('-m', '--mrna', help='file to write mRNA')
    # parser.add_argument('-c', '--cds', help='file to store CDS')
    # parser.add_argument('-p', '--protein', help='file to write translated CDS')
    parser.add_argument('-x', '--gene2mrna',
                        help='file to write gene -> mrna relations')

    args = parser.parse_args()
    MAIN(args.gxf, w=args.white, x=args.gene2mrna)
