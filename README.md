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

>Use "python align.py -h" to know more about argument types. 