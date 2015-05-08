


****************************************************************************
*******************************  RIBOABACUS  *******************************
****************************************************************************




********* Description *********


RiboAbacus is a computational tool that, given the sequence and the abundance of a transcript, provides estimation of the number of ribosomes per mRNA and transcript-specific translation efficiencies.




********* Before starting *********


Place RiboAbacus.c in the same folder of the .txt file(s) that will be used as input (see below for more details on the input files).




********* Compile and execute RiboAbacus *********


Using a Linux command-line interface: open a bash shell and type the following line to compile RiboAbacus.c and generate the executable file (that will be placed in the same folder of RiboAbacus.c)

     gcc -o RiboAbacus RiboAbacus.c -lm


then type

     ./RiboAbacus


to run RiboAbacus.




********* Settings *********


To set RiboAbacus options, an "user definition part" is present in RiboAbacus.c. The user can modified the following parameters:


CODON_USAGE 	   --->	specifies the organism codon usage bias: 1 for Homo Sapiens, 2 for Mus Musculus, 3 for Saccharomyces Cerevisiae and 0 if an input file with 				custom codon usage is given by the user (see the input files section below)
PROTEIN_INPUT 	   --->	assumes only values 0 and 1: 0 if the input file with the list of transcripts (see the input files section below) does not contain data 			related to protein amounts, 1 otherwise
RAMP_HYOUTHESIS	   --->	assumes only values 0 and 1: 0 if the ramp hypothesis is not considered, 1 otherwise
RAMP_LENGTH 	   --->	specifies the ramp length (measured in codons). If RAMP_LENGTH is set equal to 0, RAMP_LENGTH will not affect the simulation
SLOWDOWN_RATE 	   --->	specifies the ribosomes slowdown rate on the ramp (it is a percentage: for example, 70 means that the speed of the ribosomes is reduced by 				70%). If RAMP_LENGTH is set equal to 0, SLOWDOWN_RATE will not affect the simulation
WIDTH_BIN 	   --->	specifies the width of the bin for the mRNA length distribution 




********* Input files *********


The first input file must contain for each transcript TWO lines: the first reporting (in this order) the gene ID, the transcript level, the protein amount, the protein ID and the transcript ID; the second reporting the mRNA sequence. An example:

ENSG00000000971	2101.85	150.60	ENSP00000356399	ENST00000367429
ATGAGACTTCTAGCAAAGATTATTTGCCTTATGTTATGGGCTATTTGTGTAGCAGAAGATTGCAATGAACTTCCTCCAAGAAGAAATACAGAAATTCTGACAGGTTCCTGGTCTGACCAAACATATCCAGAAGGCACCCAGGCTATCTATAAATGCCGCCCTGGATATAGATCTCTTGGAAATGTAATAATGGTATGCAGGAAGGGAGAATGGGTTGCTCTTAATCCATTAAGGAAATGTCAGAAAAGGCCCTGTGGACATCCTGGAGATACTCCTTTTGGTACTTTTACCCTTACAGGAGGAAATGTGTTTGAATATGGTGTAAAAGCTGTGTATACATGTAATGAGGGGTATCAATTGCTAGGTGAGATTAATTACCGTGAATGTGACACAGATGGATGGA.......
ENSG00000002745	54.52	110.58	ENSP00000355065	ENST00000361301
ATGCAGCTCACCACTTGCCTCAGGGAGACCCTCTTCACAGGGGCTTCTCAAAAGACCTCCCTATGGTGGTTGGGCATTGCCTCCTTCGGGGTTCCAGAGAAGCTGGGCTGCGCCAATTTGCCGCTGAACAGCCGCCAGAAGGAGCTGTGCAAGAGGAAACCGTACCTGCTGCCGAGCATCCGAGAGGGCGCCCGGCTGGGCATTCAGGAGTGCGGGAGCCAGTTCAGACACGAGAGATGGAACTGCATGATCACCGCCGCCGCCACTACC.......		
ENSG00000002746	93.00	155.60	ENSP00000379228	ENST00000395891
ATGCTGCTGCACCTGTGTAGTGTGAAGAATCTGTACCAGAACAGGTTTTTAGGCCTGGCCGCCATGGCGTCTCCTTCTAGAAACTCCCAGAGCCGACGCCGGTGCAAGGAGCCGCTCCGATACAGCTACAACCCCGACCAGTTCCACAACATGGACCTCAGGGGCGGCCCCCACGATGGCGTCACCATTCCCCGCTCCACCAGCGACACTGACCTGGTCACCTCGGACAGCCGCTCCACGCTCATGGTCAGCAGCTCCTACTATTCCATCGGGCACTCTCAGGACCTGGTCATCCACTGGGACATAAAGGAGGAAGTGGACGCTGGGGACTGGATTGGCATGTACCTCATTGATGAGGTCTTGTCCGAAAACTTTCTGGACTATAAAAACCGTGGAGTCAATGGT.......

Note that if the protein amounts are not available, the corresponding column can be omitted. In this case, in the variable PROTEIN_INPUT (see the Settings section above) has to be set equal to 0.


The presence of a second input file depends on the parameter CODON_USAGE in the "user definition part" (see the Settings section above): if it is set equal to 0, an input file containing the 64 triplets and the corresponding codon usage values (without any header) must be present in the same folder of RiboAbacus.c. The following lines show an example of the organization of the file:

CUG	39.6
CCG	6.9
CAG	34.2
CGG	11.4
AUU	16
...

No specific order of the triplets is required.




********* Output files *********


RiboAbacus creates 3 output:


1)

a file in which each row corresponds to an mRNA. In each line both general transcript information and predicted data provided by RiboAbacus are reported. More in detail each line contains:
- gene ID;
- mRNA level- protein amount (optional);
- protein ID;
- transcript ID;
- transcript length;
- estimated number of ribosome per transcript;
- estimated number of ribosomes bound to the ramp region (if the ramp hypothesis was considered. Otherwise this value would be always 0);
- % of GC content of the transcript;
- % of GC content of the ramp region (if the ramp hypothesis was considered. Otherwise this value would be always 0);
- Codon Adaptation Index of the transcript;
- Codon Adaptation Index of the ramp region (if the ramp hypothesis was considered. Otherwise this value would be always 0);
- ribosome occupancy (%) based on the estimated number of ribosome per transcript;
- corrected translation efficiency based on the estimated number of ribosome per transcript;
- exit flux.


The following lines show an example of the organization of the file (ramp hypothesis: ramp length 50 codons, ribosome slowdown rate 70%):

gene ID	mRNA level	protein amount	protein ID	transcript ID	"transcript
length"	# ribosomes	"# ribosome
(ramp)"	GC%	"GC%
(ramp)"	CAI	"CAI
(ramp)"	"ribosome
occupancy"	"translation
efficiency"	exit flux
ENSG00000000971	2101.850098	150.6	ENSP00000356399	ENST00000367429	3696	32	3	39.826839	34	0.72072	0.719603	25.802345	547.0099281	44
ENSG00000004455	3493.550049	64600.76	ENSP00000346921	ENST00000354858	720	8	3	52.361111	62	0.79733	0.791008	31.911211	2287.06706	82
ENSG00000006712	8962.950195	182.54	ENSP00000221265	ENST00000221265	1596	15	3	54.699249	72	0.832251	0.887574	27.7791	2894.230964	221



2)

a file containing the frequencies of the number of ribosomes per transcript predicted by the model, useful to produce a number of ribosomes per transcript distribution that takes into account the transcript level. The following lines show an example of the organization of the file:

rib/tr	events	freq_int
0	0	0.000000
1	0	0.000000
2	1619	0.006346
3	1605	0.006291
4	2097	0.008217
5	12859	0.050386
...



3) a file containing information useful to generate a transcript length distribution with specified width of the bins and that takes into account the transcript level. The following lines show an example of the organization of the file:

length	events	freq_int
0	544	0.000710
150	1562	0.002040
300	6784	0.008862
450	8545	0.011163
600	7378	0.009639
650	8654	0.011306
...
