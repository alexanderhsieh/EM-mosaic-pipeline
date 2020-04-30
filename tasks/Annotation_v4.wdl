## Copyright Broad Institute, 2020
## This script should take as input a set of raw de novo SNVs in VCF format and annotates them using
## Ensembl Variant Effect Predictor (VEP).  This script also parses the annotation files and
## appends user-selected columns to the original raw de novo SNV VCF File.
## Output includes (1) original raw de novo SNVs annotated and (2) raw VEP output 
## NOTE: scatter-gather VEP step requires ~8GB memory (docker on mac allocates 2GB by default; need to adjust memory limit)
## 
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#Workflow Definition
###########################################################################
workflow annotate_vcf_v4 {
  
  File vcf
  String ref_ver
  String cache_dir

  File parser_script
  File parser_config
  String parser_cols  

  String final_prefix = basename(vcf, ".vcf")


  parameter_meta {
    vcf: "input raw de novo SNVs in VCF format"
    ref_ver: "reference genome version; e.g. GRCh37, GRCh38"
    cache_dir: "path to VEP cache download location"
    parser_script: "path to script used to parse and append VEP columns to original input file"
    parser_config: "path to config file for parser_script"
    parser_cols: "comma-separated string listing VEP columns to parse; default: SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
  }


  meta {
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }

  # Step *: split vcf by chromosome
  call split_vcf {
    input:
    vcf = vcf
  }


  # for each vcf chunk, run VEP to get raw annotations, then parse and append user-specific columns to each original shard vcf
  
  scatter (idx in range(length(split_vcf.out))) {
    # Step 1: generate VEP annotations
    call run_vep {
      input:
      ref = ref_ver,
      vcf = split_vcf.out[idx],
      shard = "${idx}",
      cache_dir = cache_dir
    }



    # Step 2: Parse and append VEP columns to original vcf file
    call add_vep_cols {
      input:
      original_vcf = split_vcf.out[idx],
      vep_vcf = run_vep.vep_out,
      script = parser_script,
      config = parser_config, 
      cols = parser_cols,
      shard = "${idx}"
    }


  }

  # Step 3: gather shards into final output vcf
  call gather_vep {
    input:
    vcf_shards = add_vep_cols.out,
    prefix = final_prefix
  }


  #Outputs (1) original VCF with updated annotation columns and (2) raw VEP output
  output {
    
    File cat_vcf = gather_vep.out

  }

}


###########################################################################
#Task Definitions
###########################################################################

## splits vcf by chromosome
task split_vcf {

  File vcf # input vcf
  String outprefix = basename(vcf, '.vcf')

  command <<<
    # pull header lines
    grep "^#" ${vcf} > header.txt

    # sort input vcf and bgzip
    sort -k1,1V -k2,2n ${vcf} | bgzip -c > "${vcf}.gz"

    # tabix index input vcf
    tabix -p vcf "${vcf}.gz"

    # split vcf by chromosome - use tabix -l to get all contig names from tabix index
    for i in $(tabix -l ${vcf}.gz)
    do 
      (cat header.txt; tabix ${vcf}.gz $i) > "${outprefix}.$i.vcf"
    done

    ## get full directory paths
    readlink -f *.vcf > file_full_paths.txt

  >>>

  output {
    Array[File] out = glob("*.vcf") # is the issue with the glob path?
    File filepaths = "file_full_paths.txt"
  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
  }

}


#Runs VEP in the same task to preserve environment
## NOTE: REQUIRES ~8GB memory (docker on mac allocates 2GB by default) - need to increase memory limit if running locally
## see: https://gatkforums.broadinstitute.org/wdl/discussion/11522/the-job-was-aborted-from-outside-cromwell-sporadic-failures-when-running-cromwell-locally
task run_vep {
  
  String ref # GRCh37 or GRCh38
  String cache_dir # path to location of cache files
  File vcf
  String shard

  String outprefix = basename(vcf, '.vcf')  
  String outfname = "VEP_raw.${outprefix}.${shard}.vcf"
  
  command <<<
    
    vcf_name=`basename ${vcf}`
    input_path=`dirname ${vcf}` 
    output_path=$(pwd)


    # run VEP
    ## note: may need to handle potential disagreements in VEP and cache versions via arguments
    ## used workaround described in https://gatkforums.broadinstitute.org/gatk/discussion/comment/50056#Comment_50056
    
    docker run -v $input_path:/home/vcf/ -v $output_path:/home/out/ -v ${cache_dir}:/home/vep/.vep ensemblorg/ensembl-vep:latest ./vep \
      --cache \
      --dir_cache "/home/vep/.vep" \
      --offline \
      --format vcf \
      --vcf \
      --force_overwrite \
      --assembly ${ref} \
      --input_file "/home/vcf/$vcf_name" \
      --output_file "/home/out/${outfname}" \
      --no_stats \
      --pick \
      --gencode_basic \
      --hgvs \
      --symbol \
      --transcript_version \
      --canonical \
      --biotype \
      --numbers \
      --max_af \
      --af_gnomad \
      --dir_plugins "/home/vep/.vep/Plugins/" \
      --plugin dbNSFP,/home/vep/.vep/Plugins/dbNSFP4.0a.gz,CADD_phred,MPC_score,REVEL_score

  >>>

  output {
    File vep_out = "${outfname}"
  }

}


#Parses VEP output and appends user-defined columns to the original VCF file
task add_vep_cols {
  File original_vcf
  File vep_vcf
  File script
  File config 
  String cols
  
  String shard
  String outprefix = basename(original_vcf, '.vcf')


  command <<<
    python ${script} -i ${original_vcf} -v ${vep_vcf} -c ${cols} -y ${config} -o "${outprefix}.${shard}.VEP.vcf"
  >>>

  output {
    File out = "${outprefix}.${shard}.VEP.vcf"
  }

}

#Gathers shards of VEP-annotated vcf files into a single annotated vcf
task gather_vep {

  Array[File] vcf_shards 
  String prefix

  command <<<

    while read file; do
      cat $file >> "tmp.cat.vcf"
    done < ${write_lines(vcf_shards)};

    grep "^#" "tmp.cat.vcf" > "header.txt" 

    (cat header.txt; grep -v "^#" "tmp.cat.vcf") > "${prefix}.VEP.vcf"

  >>>

  output {
    File out = "${prefix}.VEP.vcf"
  }
}

