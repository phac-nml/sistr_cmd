# 330 loci cgMLST scheme for *Salmonella*

All alleles necessary for serovar prediction are included in `cgmlst.fasta`.

## Header format

Each fasta header contains the cgMLST330 marker name and allele number delimited by `|` (pipe):

```
>{cgMLST marker name}|{allele number}
AGTC...
```

## Interpretation of results

Only alleles that match with 100% identity and query coverage are considered for cgMLST profile matching.

The closest matching cgMLST allele profile from the cgMLST allele database (`cgmlst-profiles.csv`) determines the cgMLST-derived serovar prediction.

If the input genome's cgMLST profile matches at 90% or greater similarity, then the cgMLST-derived serovar prediction will be used to refine the antigen gene-derived serovar prediction and influence the overall prediction.


