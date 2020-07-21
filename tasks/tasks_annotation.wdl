version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines tasks calling de novo SNVs from a table of single-sample gVCF paths
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
		File script
	}

	String outprefix = basename(variants, '.raw.txt')  
	String outfname = "~{outprefix}.vcf"

	command {
		set -euo pipefail

		python ~{script} -i ~{variants} -o ~{outfname}

		sort -k1,1 -k2,2n ~{outfname} | bgzip -c > "~{outfname}.gz"

		tabix -p vcf "~{outfname}.gz"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script
		String cols = "SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
	}

	String outprefix = basename(original_variants, '.raw.txt')


	command {
		python ~{script} -i ~{original_variants} -v ~{vep_vcf} -c ~{cols} -o "${outprefix}.VEP.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} "~{outprefix}.PV4.txt" 
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script
		Float cutoff_or = 3.0
		Float cutoff_p = 0.001
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript ~{script} ~{infile} ~{outprefix} ~{cutoff_or} ~{cutoff_p}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script
		Float min_vaf = 0.01
		Int size = 30000000
		Float e_fp = 0.01
		Float seq_err = 0.005
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript ~{script} ~{infile} ~{outprefix} ~{min_vaf} ~{size} ~{e_fp} ~{seq_err}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script_parse

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
		bedtools intersect -wa -wb -a "tmp.bed" -b ~{lcr} ~{map} ~{seg} -filenames -sorted > "bed.isec.out.txt"

		## parse bedtools intersect output and append relevant columns to input file
		python ~{script_parse} ~{infile} "bed.isec.out.txt" > "~{outprefix}.RR.txt"

	>>>

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
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
		File script 
		Int dist = 10# distance (in bp) used to define a "cluster"
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} ~{dist} "~{outprefix}.VC.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.VC.txt"
	}
}