# GET Aegilops speltoides data
./get_data.sh

# MERGE data sets
#AB
cat 193936_*_R1*fastq.gz >DATA_Asp_AB_R1.fastq.gz
cat 193936_*_R2*fastq.gz >DATA_Asp_AB_R2.fastq.gz
cat DATA_Asp_AB_R1.fastq.gz DATA_Asp_AB_R2.fastq.gz >DATA_Asp_AB_R1R2.fastq.gz
#0B
cat 193976_*_R1*fastq.gz >DATA_Asp_0B_R1.fastq.gz
cat 193976_*_R2*fastq.gz >DATA_Asp_0B_R2.fastq.gz
cat DATA_Asp_0B_R1.fastq.gz DATA_Asp_0B_R2.fastq.gz >DATA_Asp_0B_R1R2.fastq.gz

#CLEAN
rm 1939*fastq.gz

#---------------
#      KAT
#---------------

# COUNT
# 0B
kat sect -m 21 -N AesTR.fasta DATA_Asp_0B_R1R2.fastq.gz
# AB
kat sect -m 21 -N AesTR.fasta DATA_Asp_AB_R1R2.fastq.gz
# The k-mer counting produces the *.cvg file and a file with 
# descriptive statistics.

# COMPARE
kat comp DATA_Asp_AB_R1R2.fastq.gz DATA_Asp_0B_R1R2.fastq.gz

# PLOT
kat plot spectra-mx -o OUT_spectra kat-comp-main.mx --intersection


#---------------
# KMASKER plants
#---------------

# Constrcution of k-mer index (0B)
# These two commands construct k-mer hashes that are stored in 
# an internal repository of Kmasker plants to assists 
# reusability.
Kmasker --build --seq DATA_Asp_0B_R1R2.fastq.gz --gs 5400 --in AEGSP_tutorial_0B

# Constrcution of k-mer index (AB)
Kmasker --build --seq DATA_Asp_AB_R1R2.fastq.gz --gs 5970 --in AEGSP_tutorial_AB

# COUNT
# The following commands will perform count each k-mer of the
# input sequence and show how often it occures in the k-mer 
# hashes (kindex) of the raw data set. 
Kmasker --run --fasta AesTR.fasta --kindex AEGSP_tutorial_0B
Kmasker --run --fasta AesTR.fasta --kindex AEGSP_tutorial_AB

# COMPARE
# After running this step, look into the constructed 
# KMASKER_report_overview_statistics file. It contains a table 
# with descriptive measurements about this comparative 
# analysis.
Kmasker --run --compare --fasta AesTR.fasta --kindex AEGSP_tutorial_0B AEGSP_tutorial_AB

# PLOT
Kmasker --explore --hist --occ KMASKER_kmer_counts_KDX_AEGSP_tutorial_0B_2e2hhZDIuJ.occ --list selection.ids
Kmasker --explore --hist --occ KMASKER_kmer_counts_KDX_AEGSP_tutorial_AB_2e2hhZDIuJ.occ --list selection.ids

# Please note, that when running this tutorial on your own,
# that the process ID (here, 2e2hhZDIuJ) will change and 
# output files names may need to be adjusted in the commands.

