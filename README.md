# GPC-NatComms2021

The two custom codes were used to analyze paired-end NGS data for the paper "Scalable recombinase-based gene expression cascades".

"demultiplexing paired end reads.cpp" demultiplexes the NGS reads according to the sample identifier index and the genomic locus each read aligns to (e.g. APC, MLH1, SMAD4 and TP53).

The demultiplexed reads were used as input for CRISPRESSO2 (https://crispresso.pinellolab.partners.org/submission).

"quantify scars.cpp" allows quantification of scar sequence of the gene circuit for each stage in the cascade.
