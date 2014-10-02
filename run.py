################################################################################################################
# Alignment editing program. It removes indels and stop codons and provides an option to adjust the alignment  #
# according to user defined protein or species sequence                                                        #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################


import sys
import argparse
import textwrap
import platform
import subprocess

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Seq import _translate_str
    from Bio.Data import CodonTable
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import IUPAC
except ImportError, e:
    sys.exit("Biopython not found")

parser = argparse.ArgumentParser(prog='Alignment-Handler',
                                 version= '1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t Designed at Kimbal-Braun Lab Group, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))

group = parser.add_mutually_exclusive_group()

group.add_argument('-i', type=str, help='Enter input sequence filename')
group.add_argument('-ali', type=str, help='Enter input alignment filename')
parser.add_argument('-pkg', type=str, default = 'muscle', choices = ['muscle', 'mafft'], help='Enter the alignment package name')
parser.add_argument('-arg', type=str, default = None, help='Enter mafft argument for protein alignment')

parser.add_argument('-o', type=str, required = True, help='Enter output alignment filename')
parser.add_argument('-itype', type=str, required = True, choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'], help='Enter input alignment file format')
parser.add_argument('-otype', type=str, required = True, choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'], help='Enter output alignment file format')


parser.add_argument('-ctab', type=str, default = 'Standard', choices=['Ascidian Mitochondrial', 'SGC9', 'Coelenterate Mitochondrial', 'Protozoan Mitochondrial', 'Vertebrate Mitochondrial', 'Plant Plastid', 'Thraustochytrium Mitochondrial', 'Blepharisma Macronuclear', 'Mold Mitochondrial', 'Invertebrate Mitochondrial', 'Standard', 'Trematode Mitochondrial', 'Scenedesmus obliquus Mitochondrial', 'Euplotid Nuclear', 'Yeast Mitochondrial', 'Spiroplasma', 'Alternative Flatworm Mitochondrial', 'Ciliate Nuclear', 'SGC8', 'Alternative Yeast Nuclear', 'Hexamita Nuclear', 'SGC5', 'SGC4', 'SGC3', 'SGC2', 'SGC1', 'SGC0', 'Flatworm Mitochondrial', 'Dasycladacean Nuclear', 'Chlorophycean Mitochondrial', 'Mycoplasma', 'Bacterial', 'Echinoderm Mitochondrial'],  help='Select the codon table')


args = parser.parse_args()

if args.pkg == 'mafft' and args.arg == None:
    parser.error('-arg argument required in mafft alignment mode')

table = CodonTable.ambiguous_dna_by_id[1]


def spliter(str, num):
    '''Splits the string object'''
    return [ str[start:start+num] for start in range(0, len(str), num) ]


def translator(recordData):
    proteinSeqList = list()
    for i, rec in enumerate(recordsData):
        seqT = _translate_str(str(rec.seq), table)
        print seqT
        recordData[i].seq = Seq(seqT, IUPAC.protein)

    with open('translated.fas', 'w') as fp:
        SeqIO.write(recordsData, fp, 'fasta')


def frameCheck(records):
    for i, rec in enumerate(records):
        sequence = spliter(rec.seq, 3)
        for j, splits in enumerate(sequence):
            print splits
            if splits.count('-') == 0:
                break
            if splits.count('-') == 3:
                continue
            else:
                sequence[j] = Seq("---", generic_dna)
                break
        newSeq = Seq("", generic_dna)
        for seqData in sequence:
            newSeq = newSeq + seqData

        records[i].seq = newSeq

    return records


def alignP():
    if args.pkg == 'muscle':
        if 'Darwin' in platform.system():
            subprocess.call("./muscle/muscle -in translated.fas -out tAligned.fas", shell=True)
        else:
            subprocess.call("./muscle/muscleLinux -in translated.fas -out tAligned.fas", shell=True)

    else:
        arguments = args.arg.replace('[', '').replace(']', '')
        subprocess.call("./mafft/mafft.bat %s translated.fas > tAligned.fas" %arguments, shell=True)


def cleanAli(recordNuc):
    handleP = open('tAligned.fas', 'rU')
    records = list(SeqIO.parse(handleP, 'fasta'))
    newRecord = list()
    
    for i, rec in enumerate(records):
        nucData = [[x.id, x.seq] for x in recordNuc if x.id == rec.id][0]
        nucSeqData = spliter(nucData[1], 3)
        sequence = Seq("", generic_dna); pos = 0
        for amino in rec.seq:
            if amino == '-':
                sequence = sequence + Seq("---", generic_dna)
            else:
                sequence = sequence + nucSeqData[pos]
                pos = pos + 1

        records[i].seq = Seq(sequence, generic_dna)

    with open(args.o, 'rU') as fp:
        SeqIO.write(records, fp, args.otype)



def main():
    if args.i:
        handle = open(args.i, 'rU')
    elif args.ali:
        handle = open(args.ali, 'rU')

    records = list(SeqIO.parse(handle, args.itype))
    saveRec = records

    if args.ali:
        records = frameCheck(records)
        for i, rec in enumerate(records):
            records[i].seq = rec.seq.ungap("-")

    tSeq = translator(records)
    print records[1].seq
    alignP()
    cleanAli(records)


if __name__ == "__main__":
    main()

