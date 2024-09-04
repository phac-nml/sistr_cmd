# 1.1.3
Serovar nomenclature update after USA Cantaloupe Outbreaks in Novermber 2023. The O24 and O25 antigens would not be wet-lab typed causing the collapse of certain serovars detailed below. The selected serovar is the one that will be reported by SISTR and the other serovar in the pair will be dropped and no O24 or O25 will be reported in the antigenic formula

|Serovar pair            | Serovar selected |
|------------------------|------------------|
|Soahanina - Sundsvall   |  Sundsvall       |
|Martonos - Finkenwerder | Finkenwerder   |
|Midway - Florida      |   Florida        |
|Lindern - Charity     |   Charity        |
|Bahrenfeld - Onderstepoort | Onderstepoort |
| Schalkwijk - Moussoro  |   Schalkwijk |
|Amberg - Boecker       | Boecker      |
| Carrau - Madelia      | Carrau     |
| Chichiri - Uzaramo   | Uzaramo     |
| Poano -  Stafford    | Poano       |

### Changes of serovar assignments in `sistr/data/genomes-to-serovar.txt`

|genome accession | serovar previous | serovar current |
|-----------------|------------------|-----------------|
| SAL_DA9822AA | Soahanina |Sundsvall
| SRR1815423   | 	Soahanina | Sundsvall
| SRR2889947	| Soahanina | Sundsvall
| SRR2889992	| Soahanina | Sundsvall
| SRR3996854	|Soahanina | Sundsvall
| SRR3669910	|Soahanina | Sundsvall
| SRR3732330	|Soahanina | Sundsvall
|SRR3713652	| Soahanina | Sundsvall
| SRR3713653	|Soahanina | Sundsvall
|SRR3978444	|Soahanina | Sundsvall
| SRR2011392	|Soahanina | Sundsvall
| SRR1068363	|Soahanina | Sundsvall
|ERR161888	|Soahanina | Sundsvall
|SAL_BA5034AA	|Soahanina | Sundsvall
|SAL_EA3233AA	|Soahanina | Sundsvall
|SAL_GA9094AA	|Soahanina | Sundsvall
|SRR1158155	|Soahanina | Sundsvall
|SRR2751907	|Soahanina | Sundsvall
|SRR4237685	|Soahanina | Sundsvall
|SRR5010548	|Soahanina | Sundsvall
|09_6055	|Madelia | Carrau
|11_0879	|Madelia | Carrau
|SAL_BA1830AA	|Madelia | Carrau
|SAL_CA7979AA	|Madelia | Carrau
|SAL_DA4289AA	|Madelia | Carrau
|SAL_DA7475AA	|Madelia | Carrau
|SAL_EA4948AA	|Madelia | Carrau
|SAL_FA5821AA	|Madelia | Carrau
|SAL_HA4780AA	|Madelia | Carrau
|SAL_HA4886AA	|Madelia | Carrau
|SRR1269415	|Madelia | Carrau
|SRR1548430	|Madelia | Carrau
|SRR1805645|	Madelia | Carrau
|SRR2104612	|Madelia | Carrau
|SRR2911800	|Madelia | Carrau
|SRR3933147	|Madelia | Carrau
|SRR4098716	|Madelia | Carrau
|SRR1258654	|Madelia | Carrau
|SRR1582141	|Madelia | Carrau
|SRR4019409	|Madelia | Carrau
|SRR4244476	|Madelia | Carrau
|SRR2075023	|Madelia | Carrau
|SRR5132365	|Madelia | Carrau
|SRR5051381	|Madelia | Carrau
|SRR3743984	|Madelia | Carrau
|SRR5054238	|Madelia | Carrau
|SRR3928735	|Madelia | Carrau
|SRR1586586	|Madelia | Carrau
|SRR2976043	|Madelia | Carrau
|SRR2962333	|Madelia | Carrau
|SRR3928732	|Madelia | Carrau
|SRR3928736	|Madelia | Carrau
|SRR2962332	|Madelia | Carrau
|SAL_EA2874AA	|Bahrenfeld | Onderstepoort
|SAL_FA0525AA	|Bahrenfeld | Onderstepoort
|SRR3173783	|Bahrenfeld | Onderstepoort
|SAL_DA7014AA	|Martonos | Finkenwerder
|SRR1300569	|Martonos | Finkenwerder
|SRR1973814	|Martonos | Finkenwerder

### Changes to `Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv` antigen to serovar lookup database
Removed the following entries
1. Martonos,"6,14,24",d,"1,5",,H,FALSE,enterica
2. Midway,"6,14,24",d,"1,7",,H,FALSE,enterica
3. Lindern,"6,14,[24]",d,"e,n,x",,H,FALSE,enterica
4. Bahrenfeld,"6,14,[24]","e,h","1,5",,H,FALSE,enterica
5. Moussoro,"1,6,14,25",i,"e,n,z15",,H,FALSE,enterica
6. Amberg,"6,14,24","l,v","1,7",,H,FALSE,enterica
7. Madelia,"1,6,14,25",y,"1,7",,H,FALSE,enterica
8. Soahanina,"6,14,24",z,"e,n,x",,H,FALSE,enterica
9. Chichiri,"6,14,24","z4,z24",-,,H,TRUE,enterica
10. II 4:a:z39,"1,4,12,[27]",a,z39,,B,FALSE,salamae

# 1.1.1

* Fixed issue with sorting of BLAST results (causing cgMLST types to be different between BLAST versions). Pull request #43.

# 1.1.0

* Significant updates to SISTR antigen biomarker and cgMLST database
