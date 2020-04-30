## Copyright Broad Institute, 2020
## This script should take as input a set of filtered and QC'd de novo SNVs and return (1) the same de novo file
## with additional columns for p-value, likelihood ratio, and posterior odds (2) a separate file containing
## only candidate mosaic variants (variants with posterior odds > user-supplied cutoff) and (3) set of plots 
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
workflow EM_mosaic {
  
  File dnsnvs 
  File em_mosaic
  Int postcut
  
  parameter_meta{
    dnsnvs: "de novo SNVs file with cols {id, chr, pos, ref, alt, refdp, altdp} at minimum"
    em_mosaic: "full path to EM-mosaic R script"
    postcut: "posterior odds cutoff score (default: 10) for defining mosaic vs. germline"
  }
  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }

  #run EM_mosaic
  call detect_mosaic {
  	input:
  	infile = dnsnvs,
    script = em_mosaic,
    postcut = postcut
  }

  #Outputs (1) de novos with score columns (2) candidate mosaics (3) 5x plots
  output {
      File denovos = detect_mosaic.denovos
      File mosaics = detect_mosaic.mosaics
      File plot_EM = detect_mosaic.plot_EM
      File plot_QQ = detect_mosaic.plot_QQ
      File plot_overdispersion = detect_mosaic.plot_overdispersion
      File plot_dp_vs_vaf = detect_mosaic.plot_dp_vs_vaf
      File plot_vaf_vs_post = detect_mosaic.plot_vaf_vs_post
    }

}


###########################################################################
#Task Definitions
###########################################################################

#Scores variants in dnSNVs file
task detect_mosaic {
  File infile
  File script
  Int postcut 
  String outprefix = basename(infile, '.txt')

  command {
    Rscript ${script} ${infile} ${outprefix} ${postcut} 
  }

  output {
    File denovos = "${outprefix}.denovo.txt"
    File mosaics = "${outprefix}.candidates.txt"
    File plot_EM = "${outprefix}.EM.pdf"
    File plot_QQ = "${outprefix}.QQ.pdf"
    File plot_overdispersion = "${outprefix}.overdispersion.pdf"
    File plot_dp_vs_vaf = "${outprefix}.dp_vs_vaf.pdf"
    File plot_vaf_vs_post = "${outprefix}.vaf_vs_post.pdf"
  }
}