# mod_seq
This package requires installation of the ShapeMapper package from Kevin Weeks' lab.
Download ShapeMapper from http://www.chem.unc.edu/rna/software.html compile it, and add it to your PATH.

note: increasing maxProc in shapemapper.py in shapemapper package improves parrallelization when run locally with multiple
cores available.

cutadapt is required: http://cutadapt.readthedocs.org/en/stable/guide.html

fastX toolkit is required: http://hannonlab.cshl.edu/fastx_toolkit/

the bokeh interactive plotting package is optionally required:
http://bokeh.pydata.org/en/latest/docs/user_guide/quickstart.html#userguide-quickstart


explanation of settings files:
since I can't figure out how to add comments to the settings files, please see the explanation below, as well as the example file.


[ input ]
fastq_dir = /Users/boris/Green_Lab/Book_2/2.47/E2.45_sequencing/FASTQ #full path to folder with gzipped fastq files
#file names for fastq.gz files to be analyzed, within above folder
fastq_gz_files = ["HJLLNBCXX_1_AGGTTT_1.fastq.gz","HJLLNBCXX_1_CCTGAG_1.fastq.gz","HJLLNBCXX_1_GAACCC_1.fastq.gz"]
#sample names for each fastq.gz file, in same order
sample_names = ["80S_no_DMS","80S","80S_DMSO"]
#names of what are considered "experimental" samples, which will be normalized to a control
experimentals = ["80S_DMSO",]
#names of unmodified samples, that are not treated with DMS or other agent, in matching order to above
no_mod_controls= ["80S_no_DMS""]
#names of modified control samples, that are treated with DMS or other agent, but not the experimental condition, in matching order to above
with_mod_controls  = ["80S"]

[ parameters ]
experiment_name = E2.47
#a FASTA file with the RNAs (assumed ribosomal, but not required) to map to
rrna_fasta = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/rRNA_Sequences_ucsc_20110829.txt
#the base file that will be modified by the pipeline to configure a shapemapper run. the default should usually work
shapemapper_ref_file = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/shapemapper_BASE.cfg
#for making ROC curves, true positive and true negative calls for 18S and 25S rRNA
tptn_file_18s = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/structure_ROC_curves/18S_4V88_2.0_AC.txt
tptn_file_25s = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/structure_ROC_curves/25S_4V88_2.0_AC.txt

#nucleotides to exclude that will always show up as modified, even without modifying agent
exclude_constitutive = {"S.c.18S_rRNA":[1191],"S.c.25S__rRNA":[645,2142,2634,2843]}
#1-p-value cutoff. Include the mutliple testing correction here
confidence_interval_cutoff = 0.999999
#fold change required, on top of p-value pass, to be considered protected or deprotected
fold_change_cutoff = 2

#base scripts used to create annotated pymol scripts
pymol_base_script = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/structure_highlighting/ban_nomenclature_S_cerevisiae_80S_PyMOL_highlighting_base.txt
pymol_base_script_colorchange = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/structure_highlighting/ban_nomenclature_S_cerevisiae_80S_PyMOL_colored_by_change_base.txt

#tsv of groups of nucleotides with relevant functions (like subunit interface, or A site)
functional_groupings = /Users/boris/Green_Lab/Book_2/2.36/mod_seq/structure_highlighting/S.cerevisiae_rRNA_functional_groups.txt

#setting this to Ture requires Bokeh package
make_interactive_plots = True
#Whatever was ligated to RNA 3' end in library prep
adaptor_sequence = CACTCGGGCACCAAGGAC
#what nucleotides your probe hits
affected_nucleotides = AC
trim_adaptor = True
discard_untrimmed = True
min_post_adaptor_length = 40
first_base_to_keep = 4
last_base_to_keep = 100
min_base_quality = 20
min_mapping_quality = 10
#for troubleshooting
collapse_identical_reads = False
force_read_resplit = False
force_recollapse = False
force_retrim = False
force_remapping = False
force_index_rebuild = False
force_shapemapper = False
force_recount = True
#rarely used
error_dir = /Users/boris/Green_Lab/Book_2/2.47/analysis/error/

[ output ]
#where to put the results
results_dir = /Users/boris/Green_Lab/Book_2/2.47/analysis
