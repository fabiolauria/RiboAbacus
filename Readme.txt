


*******************************  RIBOABACUS  *******************************





********* Description *********


RiboAbacus is a computational tool that, given the sequence and the abundance of a transcript, provides estimation of the number of ribosomes per mRNA and transcript-specific translation efficiencies.




********* Before starting *********


Place RiboAbacus.c in the same folder of the .txt file that will be used as input (see below for more details on the input file).




********* Compile and execute RiboAbacus *********


Using a Linux command-line interface: open a command shell window and type the following line to obtain the executable file

> gcc -o RiboAbacus RiboAbacus.c -lm


then type

./RiboAbacus


to compile RiboAbacus.




********* Input file *********


The input file should contain for each transcript TWO lines: the first reporting (in this order) the gene ID, the transcript level, the protein amount, the protein ID and the transcript ID; the second reporting the mRNA sequence. An example:

ENSG00000000971	2101.85	150.60	ENSP00000356399	ENST00000367429
ATGAGACTTCTAGCAAAGATTATTTGCCTTATGTTATGGGCTATTTGTGTAGCAGAAGATTGCAATGAACTTCCTCCAAGAAGAAATACAGAAATTCTGACAGGTTCCTGGTCTGACCAAACATATCCAGAAGGCACCCAGGCTATCTATAAATGCCGCCCTGGATATAGATCTCTTGGAAATGTAATAATGGTATGCAGGAAGGGAGAATGGGTTGCTCTTAATCCATTAAGGAAATGTCAGAAAAGGCCCTGTGGACATCCTGGAGATACTCCTTTTGGTACTTTTACCCTTACAGGAGGAAATGTGTTTGAATATGGTGTAAAAGCTGTGTATACATGTAATGAGGGGTATCAATTGCTAGGTGAGATTAATTACCGTGAATGTGACACAGATGGATGGA.......
ENSG00000002745	54.52	110.58	ENSP00000355065	ENST00000361301
ATGCAGCTCACCACTTGCCTCAGGGAGACCCTCTTCACAGGGGCTTCTCAAAAGACCTCCCTATGGTGGTTGGGCATTGCCTCCTTCGGGGTTCCAGAGAAGCTGGGCTGCGCCAATTTGCCGCTGAACAGCCGCCAGAAGGAGCTGTGCAAGAGGAAACCGTACCTGCTGCCGAGCATCCGAGAGGGCGCCCGGCTGGGCATTCAGGAGTGCGGGAGCCAGTTCAGACACGAGAGATGGAACTGCATGATCACCGCCGCCGCCACTACC.......		
ENSG00000002746	93.00	155.60	ENSP00000379228	ENST00000395891
ATGCTGCTGCACCTGTGTAGTGTGAAGAATCTGTACCAGAACAGGTTTTTAGGCCTGGCCGCCATGGCGTCTCCTTCTAGAAACTCCCAGAGCCGACGCCGGTGCAAGGAGCCGCTCCGATACAGCTACAACCCCGACCAGTTCCACAACATGGACCTCAGGGGCGGCCCCCACGATGGCGTCACCATTCCCCGCTCCACCAGCGACACTGACCTGGTCACCTCGGACAGCCGCTCCACGCTCATGGTCAGCAGCTCCTACTATTCCATCGGGCACTCTCAGGACCTGGTCATCCACTGGGACATAAAGGAGGAAGTGGACGCTGGGGACTGGATTGGCATGTACCTCATTGATGAGGTCTTGTCCGAAAACTTTCTGGACTATAAAAACCGTGGAGTCAATGGT.......



Note that if the protein amounts are not available the entire corresponding column could be omitted. In this case, in the variable PROTEIN_INPUT (see the "user definition part" in RiboAbacus.c) has to be set equal to 0.




********* Output file *********


RiboAbacus will return a .txt file in which each row corresponds to an mRNA. In each line both general transcript information and predicted data provided by RiboAbacus are reported. More in detail each line contains:
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


An example (ramp hypothesis: ramp length 50 codons, ribosome slowdown rate 70%):

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


























