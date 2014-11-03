#Alignment-Handler
Performs codon alignment by using muscle/mafft package.
It automatically identifies a correct frame to perform codon alignment.

##Instructions
Use 
```
python align.py -i inputFileName -itype inputAlignMentType -o OutPutFileName -otype outputAlignMentType
```

You can choose to ignore stop codon check using -ign argument.

```
python align.py -i inputFileName -itype inputAlignMentType -o OutPutFileName -otype outputAlignMentType -ign
```

####Handling low quality CDS 
Several coding sequences contains in-frame stop codons and are considered as low quality sequences. 
Alignment-Handler '-omit' argument screens these sequences and omits it from the alignment dataset.

```
python align.py -i inputFileName -itype inputAlignMentType -o OutPutFileName -otype outputAlignMentType -omit
```

>Use "python align.py -h" to know more about argument types. 