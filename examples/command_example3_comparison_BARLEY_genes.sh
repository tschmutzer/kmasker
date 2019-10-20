
# GET DATA
wget -c https://webblast.ipk-gatersleben.de/barley_ibsc/downloads/160517_Hv_IBSC_PGSB_r1_CDS_HighConf_REPR_annotation.fasta.gz
gunzip 160517_Hv_IBSC_PGSB_r1_CDS_HighConf_REPR_annotation.fasta.gz

# MAIN ANALYSIS
# First, we call the main command using pre-constructed k-mer hash structures 
# from two WGS data sets (Morex and Igri)
Kmasker --run --compare --fasta 160517_Hv_IBSC_PGSB_r1_CDS_HighConf_REPR_annotation.fasta --kindex BARLEY_morex_tutorial BARLEY_igri_tutorial --strict

# DIFFERENT
# In total, 49 genes have at least one sequence segment of length larger 
# than 30 bp, that has significantly increased k-mer values in morex compared
# to Igri.
grep KDC KMASKER_diverse_regions_KDX_BARLEY_igri_tutorial_4uAG4sImvL.gff|grep -v KDR|awk -F"\t" '{if(($5-$4)>=30) print $0}'

# ABSENT
# In total, 232 genes are found with large parts absent in winter barley Igri. 
awk -F"\t" '{if(($16>=80)&&($14<20)) print $0}' KMASKER_report_statistics_compare_*.txt >Gene_candidates.txt

# PLOT
# In compare mode Kmasker plants generates a report that contains several 
# descriptive measurements and which can be visualized to compare the two 
# datasets.
cut -f 1,14,16 KMASKER_report_statistics_compare_* >INPUT_plotting.txt
Kmasker --explore --hexplot --file INPUT_plotting.txt

# To get an overview about the general distribution the violin plot can assist.
Kmasker --explore --violin  --file INPUT_plotting.txt
# The violin plot shows at the x-axis the percentage of the gene length that is
# missing. The tendency is that Morex (right side) has more complete genes compared
# to Igri (left side). Also, the peak at 100, which represents those genes where
# 100% of the gene length is absent, is indicative. Here, based on k-mer counts 
# Igri (31) has more absent genes than Morex (15).
