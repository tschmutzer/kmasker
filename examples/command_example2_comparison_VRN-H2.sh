#--------------
# GET DATA
---------------

./get_igri_PE.sh
./get_morex_PE.sh

#---------------
#      KAT
#---------------

# COUNT
# IGRI
kat sect -m 21 -N -o KAT_Igri sequence_and_primer_VRN-H2.fasta R1R2_Sample_Igri.fastq.gz
# MOREX
kat sect -m 21 -N -o KAT_Morex sequence_and_primer_VRN-H2.fasta R1R2_Sample_Morex.fastq.gz

# COMPARE
kat comp -o KAT_compare R1R2_Sample_Igri.fastq.gz R1R2_Sample_Morex.fastq.gz

# PLOT
kat plot spectra-mx -o OUT_spectra KAT_compare_kat-comp-main.mx --intersection


#---------------
# KMASKER plants
#---------------

# Construction of k-mer index (IGRI)
Kmasker --build --seq R1_Sample_Igri.fastq.gz R2_Sample_Igri.fastq.gz --gs 5100 --in BARLEY_igri_tutorial
# Construction of k-mer index (IGRI)
Kmasker --build --seq R1_Sample_Morex.fastq.gz R2_Sample_Morex.fastq.gz --gs 5100 --in BARLEY_morex_tutorial

# COUNT
Kmasker --run --fasta sequence_and_primer_VRN-H2.fasta --kindex BARLEY_igri_tutorial
Kmasker --run --fasta sequence_and_primer_VRN-H2.fasta --kindex BARLEY_morex_tutorial

# COMPARE
Kmasker --run --compare --fasta sequence_and_primer_VRN-H2.fasta --kindex BARLEY_igri_tutorial BARLEY_morex_tutorial
