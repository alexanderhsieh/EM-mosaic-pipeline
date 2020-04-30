## Copyright Broad Institute, 2020
## This script takes as input a set of raw de novo SNV calls and runs them through a filtering pipeline
## Filter codes include:
##  SB = separate test for strand bias using adfref, adfalt, adrref, adralt values
##  FDR = FDR-based minimum altdp filter
##  RR = repetitive regions (Mappability < 1, Segmental Duplication > 0.95, Low-Complexity Regions)
##  VC = variant clusters (two variants within x bp of each other)
##  MAF = population minor allele frequency (gnomAD)
##  COD = protein-coding
##  MUC-HLA = MUC- and HLA- genes; tend to be very noisy
##  CAF = cohort allele frequency
##  IGV = IGV manual review result
##  dbSNP = fail if variant present in dbSNP
##  OUT = outlier_flag result
##  
## Variants passing all filters are then passed to EM-mosaic for mosaic variant discovery
##  
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#WORKFLOW DEFINITION
###########################################################################
workflow ApplyFilters {

  File anno 

  File script_pv4 

  File script_sb 
  Float sb_or 
  Float sb_p 

  File script_fdr 
  Int fdr_min 
  Int fdr_size 
  Float fdr_e_fp 
  Float fdr_seq_err 

  File script_rr_bash 
  File script_rr_parse 
  String rr_bedtools_binpath 
  File rr_map 
  File rr_seg 
  File rr_lcr 

  File script_vc 
  Int vc_dist 

  File script_update_filter_col
  File update_filter_col_config 

  File script_outlier 
  Int cohort_size
  Int cutoff

  File script_filtct 

  File script_printpass 

  parameter_meta {
    anno: "full path to file containing annotated denovos from ANNOTATION step"
    script_pv4: "full path to script filter_PV4.py"
    script_sb: "full path to script filter_strandbias.R"
    sb_or: "strand bias test Fisher's Exact Test OR threshold value"
    sb_p: "strand bias test Fisher's Exact Test p-value threshold"
    script_fdr: "full path to script filter_fdrmin.R"
    fdr_min: "floor value for FDR-based minimum Nalt threshold calculation"
    fdr_size: "region size for FDR-based min Nalt threshold calculation"
    fdr_e_fp: "threshold for expected number of false positives across entire region size"
    fdr_seq_err: "sequencing error rate for FDR-based minimum Nalt threshold calculation"
    script_rr_bash: "full path to bash script run_bedtools.sh"
    script_rr_parse: "full path to python script parse_bedtools_isec.py"
    rr_bedtools_binpath: "full path to folder containing bedtools executable"
    rr_map: "full path to mappability BEDfile"
    rr_seg: "full path to segmental duplication BEDfile"
    rr_lcr: "full path to LCR BEDfile"
    script_vc: "full path to script flag_vclust.py"
    vc_dist: "distance in bp to define variant clusters"
    script_update_filter_col: "full path to script update_filter_col.py"
    update_filter_col_config: "full path to config .yml file for update_filter_col.py"
    script_outlier: "full path to script filter_outlier.R"
    cohort_size: "number of samples in cohort"
    cutoff: "default=1; use the point where the number of samples with a given number of de novos falls below cutoff"
    script_filtct: "full path to script get_filt_ct.py"
    script_printpass: "full path to script print_pass_only.py"
  }
  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }


  #run PV4 filter
  call filter_PV4 {
  	input:
    infile = anno,
  	script = script_pv4
  }
  #run SB (strand bias) filter
  call filter_SB {
    input:
    infile = filter_PV4.outfile,
    script = script_sb,
    cutoff_or = sb_or,
    cutoff_p = sb_p
  }
  #run FDR (FDR-based min altdp) filter
  call filter_FDR {
    input:
    infile = filter_SB.outfile,
    script = script_fdr,
    min = fdr_min,
    size = fdr_size,
    e_fp = fdr_e_fp,
    seq_err = fdr_seq_err
  }
  #run RR (repeat region) filter
  call filter_RR {
    input:
    infile = filter_FDR.outfile,
    script_bash = script_rr_bash,
    script_parse = script_rr_parse,
    bedpath = rr_bedtools_binpath,
    map = rr_map,
    seg = rr_seg,
    lcr = rr_lcr
  }
  #run VC (variant cluster) filter
  call filter_VC {
    input:
    infile = filter_RR.outfile,
    script = script_vc,
    dist = vc_dist
  }
  #run update_filter_column script to combine filter flags into single column
  call update_filt_col {
    input:
    infile = filter_VC.outfile,
    script = script_update_filter_col,
    config = update_filter_col_config
  }
  #run OUT (outlier) filter
  call filter_OUT{
    input:
    infile = update_filt_col.outfile,
    script = script_outlier,
    cohort_size = cohort_size,
    cutoff = cutoff
  }
  #run script to summarize counts of variants flagged by each filter
  call summarize_counts{
    input:
    infile = filter_OUT.outfile,
    script = script_filtct
  }
  #run script to write out variants passing all filters, to be used as input to EM-mosaic
  call print_pass_vars{
    input:
    infile = filter_OUT.outfile,
    script = script_printpass
  }


  #Outputs 
  output {
      File denovos_all = filter_OUT.outfile # full denovos table with all annotations
      File denovos_PASS = print_pass_vars.outfile # only denovos passing all filters, to be passed to EM-mosaic scoring step
      File filter_summary = summarize_counts.outfile # summary table detailing how many variants failed each filter step
    }

}





###########################################################################
#Task Definitions
###########################################################################

#Adds PV4 columns to variants file for downstream filtering
task filter_PV4 {
  File infile
  File script
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} "${outprefix}.PV4.txt" 
  }

  output {
    File outfile = "${outprefix}.PV4.txt"
  }
}

#Adds strand_bias columns to variants file for downstream filtering
task filter_SB {
  File infile
  String outprefix = basename(infile, '.txt')
  File script
  Float cutoff_or
  Float cutoff_p

  command {
    Rscript ${script} ${infile} ${outprefix} ${cutoff_or} ${cutoff_p}
  }

  output {
    File outfile = "${outprefix}.SB.txt"
  }
}

#Adds FDR-based minimum Nalt filter-related columns to variants file for downstream filtering
task filter_FDR {
  File infile
  String outprefix = basename(infile, '.txt')
  File script
  Int min 
  Int size 
  Float e_fp
  Float seq_err 

  command {
    Rscript ${script} ${infile} ${outprefix} ${min} ${size} ${e_fp} ${seq_err}
  }

  output {
    File outfile = "${outprefix}.FDR.txt"
  }
}

#Adds Repeat Region (Mappability, SegDup, LCR) filter-related columns to variants file for downstream filtering
task filter_RR {
  File infile
  File script_bash
  File script_parse
  String bedpath
  File map
  File seg
  File lcr
  String outprefix = basename(infile, '.txt')

  command {
    sh ${script_bash} -v ${infile} -b ${bedpath} -m ${map} -s ${seg} -l ${lcr} -p ${script_parse}
  }

  output {
    File tmpbed = "tmp.bed" # temporary BEDfile created from variants file
    File isec_file = "bed.isec.out.txt" # BEDtools intersect output
    File outfile = "${outprefix}.RR.txt"
  }
}

#Adds Variant Cluster filter-related columns to variants file for downstream filtering
task filter_VC {
  File infile
  File script 
  Int dist # distance (in bp) used to define a "cluster"
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} ${dist} "${outprefix}.VC.txt"
  }

  output {
    File outfile = "${outprefix}.VC.txt"
  }
}

#Parses filter information from upstream steps and combines into an updated filter column
#Note: assumes 'filter' col already exists (from upstream reformat_vcf.py step)
#Note: also applies filters for 
#   (1) popfreq (MAF)
#   (2) protein-coding (COD)
#   (3) MUC/HLA genes (MUC-HLA)
#   (4) dbSNP (dbSNP)
task update_filt_col {
  File infile 
  File script 
  File config 

  command {
    python ${script} -i ${infile} -c ${config} -o "ADfile.filt.txt"
  }

  output {
    File outfile = "ADfile.FILT.txt"
  }
}

#Adds Outlier filter-related columns to variants file for downstream filtering
task filter_OUT {
  File infile
  File script
  Int cohort_size
  Int cutoff
  String outprefix = basename(infile, '.txt')

  command {
    Rscript ${script} ${infile} ${outprefix} ${cohort_size} ${cutoff}
  }

  output {
    File outfile = "${outprefix}.OUT.txt"
  }
}

#Parses filter column and summarize how many variants are flagged by each filter
task summarize_counts {
  File infile 
  File script 
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} "${outprefix}.SUMMARY_COUNTS.txt"
  }

  output {
    File outfile = "${outprefix}.SUMMARY_COUNTS.txt"
  }
}

#Prints out all variants that pass all filters and that belong to samples that are not flagged as outliers
task print_pass_vars {
  File infile 
  File script 
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} "${outprefix}.PASS.txt"
  }

  output {
    File outfile = "${outprefix}.PASS.txt"
  }
}
