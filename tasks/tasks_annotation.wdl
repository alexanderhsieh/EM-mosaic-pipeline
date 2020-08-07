version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines tasks for annotating an raw de novo SNV callset
##
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# converts tab-separated variants file into vcf format for VEP
task txt_to_vcf {
	input {
		File variants
	}

	String outprefix = basename(variants, '.raw.txt')  
	String outfname = "~{outprefix}.vcf"

	command {
		set -euo pipefail

		python /opt/convert_txt_to_vcf.py -i ~{variants} -o ~{outfname}

		sort -k1,1 -k2,2n ~{outfname} | bgzip -c > "~{outfname}.gz"

		tabix -p vcf "~{outfname}.gz"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outfname}.gz"
		File idx = "~{outfname}.gz.tbi"
	}

}


#Runs VEP in the same task to preserve environment
## NOTE: REQUIRES ~8GB memory (docker on mac allocates 2GB by default) - need to increase memory limit if running locally
## see: https://gatkforums.broadinstitute.org/wdl/discussion/11522/the-job-was-aborted-from-outside-cromwell-sporadic-failures-when-running-cromwell-locally
task run_vep {
	input {
		String ref # GRCh37 or GRCh38
		File cache_dir # path to location of cache files
		String cache_version = "100"
		File vcf
		File vcf_idx
		Int disk_size = 125 # test 100G for VEP?
	}


	String cache_dirname = basename(cache_dir, '.tar.gz')
	#String outprefix = basename(vcf, '.vcf')  
	String outprefix = basename(vcf, '.vcf.gz')  
	String outfname = "VEP_raw.~{outprefix}.vcf"

	command {

		echo "## extracting and localizing cache directory"
		time tar -xf ~{cache_dir} -C ~/
		CACHE_PATH=$(cd ~/~{cache_dirname}; pwd)
		echo "## success; see cache directory -- "$CACHE_PATH

		/opt/vep/src/ensembl-vep/vep \
		--cache \
		--dir_cache "$CACHE_PATH" \
		--cache_version ~{cache_version} \
		--offline \
		--fork 4 \
		--format vcf \
		--vcf \
		--assembly ~{ref} \
		--input_file ~{vcf} \
		--output_file ~{outfname} \
		--no_stats \
		--pick \
		--symbol \
		--canonical \
		--biotype \
		--max_af 


		rm -rf ~{cache_dir}
		rm -rf $CACHE_PATH

	}

	runtime {
		docker: "ensemblorg/ensembl-vep:latest"
		disks: "local-disk " + disk_size + " HDD"
		bootDiskSizeGb: disk_size
		memory: "16G"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File vep_out = "~{outfname}"
	}

}

#Parses VEP output and appends user-defined columns to the original VCF file
task add_vep_cols {
	input {
		File original_variants
		File vep_vcf
		String cols = "SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
	}

	String outprefix = basename(original_variants, '.raw.txt')


	command {
		python /opt/parse_and_append_vep_cols.py -i ~{original_variants} -v ~{vep_vcf} -c ~{cols} -o "${outprefix}.VEP.txt"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.VEP.txt"
	}

}

#Adds PV4 columns to variants file for downstream filtering
task flag_PV4 {
	input {
		File infile
	}

	String outprefix = basename(infile, '.txt')

	command {
		python /opt/filter_GATK_RankSum.py ~{infile} "~{outprefix}.PV4.txt" 
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.PV4.txt"
	}
}

#Adds strand_bias columns to variants file for downstream filtering
task flag_SB {
	input {
		File infile
		Float cutoff_or = 3.0
		Float cutoff_p = 0.001
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript /opt/filter_strandbias.R ~{infile} ~{outprefix} ~{cutoff_or} ~{cutoff_p}
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.SB.txt"
	}
}

#Adds FDR-based minimum Nalt filter-related columns to variants file for downstream filtering
task flag_FDR {
	input {
		File infile
		Float min_vaf = 0.01
		Int size = 30000000
		Float e_fp = 0.01
		Float seq_err = 0.005
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript /opt/filter_fdrmin.R ~{infile} ~{outprefix} ~{min_vaf} ~{size} ~{e_fp} ~{seq_err}
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.FDR.txt"
	}
}

#Adds Repeat Region (Mappability, SegDup, LCR) filter-related columns to variants file for downstream filtering
task flag_RR {
	input {
		File infile

		File map
		File seg
		File lcr

		Int disk_size = 100
	}

	String outprefix = basename(infile, '.txt')

	command <<<

		set -euo pipefail

		BEDPATH=`which bedtools`
		echo $BEDPATH > 'bedpath.txt'

		## NOTE: using alternate command brackets to accommodate awk - otherwise will get unrecognized token error
		## format input as bedfile
		awk -F '\t' '{if($1!="id") print "chr"$2"\t"$3"\t"$3}' ~{infile} | sort -k1,1 -k2,2n > "tmp.bed"

		## run bedtools intersect
		bedtools intersect -wa -wb -a "tmp.bed" -b ~{lcr} ~{map} ~{seg} -filenames > "bed.isec.out.txt"

		## parse bedtools intersect output and append relevant columns to input file
		python /opt/parse_bedtools_isec.py ~{infile} "bed.isec.out.txt" > "~{outprefix}.RR.txt"

	>>>

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
		memory: "16G"
		disks: "local-disk " + disk_size + " HDD"
		bootDiskSizeGb: disk_size

	}

	output {
		File bedpath = "bedpath.txt"
		File tmpbed = "tmp.bed" # temporary BEDfile created from variants file
		File isec_file = "bed.isec.out.txt" # BEDtools intersect output
		File out = "~{outprefix}.RR.txt"
	}
}


#Adds Variant Cluster filter-related columns to variants file for downstream filtering
task flag_VC {
	input {
		File infile
		Int dist = 10# distance (in bp) used to define a "cluster"
	}

	String outprefix = basename(infile, '.txt')

	command {
		python /opt/flag_vclust.py ~{infile} ~{dist} "~{outprefix}.VC.txt"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.VC.txt"
	}
}

# estimates cohort AF
task estimate_cohort_AF{

	input {
		File infile
		Int cohort_size
	}

	String outprefix = basename(infile, '.denovo.VEP.PV4.SB.FDR.RR.VC.txt')
	String outfname = "AC.~{outprefix}.txt"

	command {
		python /opt/estimate_cohort_AF-SIMPLE.py -i ~{infile} -n ~{cohort_size} -o ~{outfname}
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outfname}"
	}
}


#Adds cohort allele frequency filter-related columns to variants file for downstream filtering
task flag_CAF {

	input {
		File infile
		File caf_file
	}

	String outprefix = basename(infile, '.txt')

	command {
		python /opt/filter_cohort_AF.py -i ~{infile} -c ~{caf_file} -o "~{outprefix}.CAF.txt"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.CAF.txt"
	}

}

#Adds Outlier filter-related columns to variants file for downstream filtering
task flag_outlier {

	input {
		File infile
		Int cohort_size
		Int exp
		Int cutoff
	}

	String outprefix = basename(infile, '.txt')

	command {
		Rscript /opt/filter_outlier_v2.R ~{infile} ~{outprefix} ~{cohort_size} ~{exp} ~{cutoff}
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.OUT.txt"
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

	input {
		File infile 
	}

	String outsuffix = basename(infile, '.VEP.PV4.SB.FDR.RR.VC.CAF.OUT.txt')
	String outfname = "ADfile.~{outsuffix}.txt"

	command {
		python /opt/update_filter_col_v2.py -i ~{infile} -o ~{outfname}
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outfname}"
	}

}


#Parses filter column and summarize how many variants are flagged by each filter
task summarize_counts {

	input {
		File infile 
	}

	String outprefix = basename(infile, '.txt')

	command {
		python /opt/get_filt_ct_v2.py ~{infile} "~{outprefix}.SUMMARY_COUNTS.txt"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outprefix}.SUMMARY_COUNTS.txt"
	}
}

#Prints out all variants that pass all filters and that belong to samples that are not flagged as outliers
task print_pass_vars {

	input {
		File infile 
	}
	
	String outprefix = basename(infile, '.txt')

	command {
		python /opt/print_pass_only.py ~{infile} "~{outprefix}.PASS.txt"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:7Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outprefix}.PASS.txt"
	}
}