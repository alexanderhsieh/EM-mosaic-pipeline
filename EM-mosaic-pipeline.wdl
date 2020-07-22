version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines the workflow for using EM-mosaic to call mosaic SNVs from 
## a cohort of trios (in single-sample gzipped gVCF format).  
## Note: currently limited to human whole-exome-sequencing data
## 
## Requirements:
## - cohort of trios with gVCFs for each individual
## - pedigree file with columns: family ID, sample ID, father ID, mother ID, sex, affected status
## - sample map (from Picard) with 2 columns: sample ID, GVCF path
## - sample table (from Terra) with columns: entity:sample_id, output_vcf, output_vcf_index
## - 
## 
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_gvcf_to_denovo.wdl" as gvcf_to_denovo
import "https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_annotation.wdl" as annotation
import "https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_filtered_denovo_to_mosaic.wdl" as call_mosaic


###########################################################################
# WORKFLOW DEFINITION
###########################################################################
workflow EM_mosaic_pipeline {

	input {

		## gvcf to denovo
		File sample_table
		File sample_map
		File ped
		File ref_fasta
		File ref_fasta_index
		Float pb_min_vaf
		Int par_min_dp
		Int par_max_alt
		String output_prefix
		String output_suffix

		## annotation
		String ref_ver
		File cache_dir
		File rr_map 
		File rr_seg 
		File rr_lcr 

		## apply filters
		String output_prefix
		String CAF_outprefix
		Int cohort_size
		Int expected_dnsnvs 
		Int case_cutoff 

		## filtered denovo to mosaic
		Int postcut
		Int cohort_size 
		String output_prefix

	}

	##########################################################################
	## Call de novo SNVs from gVCFs in sample map
	##########################################################################

	## Parse sample table downloaded from Terra (.tsv)
	## Generate table that enables grouping trio gvcfs together as Array[File]
	call gvcf_to_denovo.read_terra_table{
		input:
			terra_table = sample_table,
			sample_map = sample_map,
			ped = ped
	}

	## Read in table and coerce as Array[File] corresponding to trios
	## output 7-column table: 
	## 		sample_id, 
	##      sample gvcf google bucket path, father gvcf path, mother gvcf path
	##      sample gvcf index google bucket path, father gvcf index path, mother gvcf index path
	call gvcf_to_denovo.read_table {
		input:
			table = read_terra_table.out
	}

	## i index (rows) correspond to individual samples
	scatter (i in range(length(read_table.out))) {

		Int n_cols = length(read_table.out[0])

		## parse columns containing gvcf google bucket paths 
		scatter (j in range(n_cols)) {
			if (j >=1 && j<=3) {
				File gvcf_columns = read_table.out[i][j]
			}
		}

		## parse columns containing gvcf index google bucket paths
		scatter (j in range(n_cols)) {
			if (j >=4 && j<=6) {
				File gvcf_index_columns = read_table.out[i][j]
			}
		}

		String sample_id = read_table.out[i][0]
		Array[File] selected_gvcf_columns = select_all(gvcf_columns)
		Array[File] selected_gvcf_index_columns = select_all(gvcf_index_columns)

		call gvcf_to_denovo.merge_trio_gvcf {
			input:
				sample_id = sample_id,
				trio_gvcf_array = selected_gvcf_columns,
				trio_gvcf_index_array = selected_gvcf_index_columns,
				ref_fasta = ref_fasta,
				ref_fasta_index = ref_fasta_index
		}

		call gvcf_to_denovo.call_denovos {
			input:
				sample_id = sample_id,
				trio_readgroup_ids = merge_trio_gvcf.rg_ids,
				gvcf = merge_trio_gvcf.out_gvcf,
				pb_min_vaf = pb_min_vaf,
				par_max_alt = par_max_alt,
				par_min_dp = par_min_dp,
				output_suffix = output_suffix
		}

	}

	call gvcf_to_denovo.gather_shards {
		input:
			shards = call_denovos.out,
			headers = call_denovos.head,
			output_prefix = output_prefix,
			output_suffix = output_suffix
	}

	##########################################################################
	## Annotate raw de novo variants
	##########################################################################
	
	# convert txt to vcf
	call annotation.txt_to_vcf {
		input:
			variants = gather_shards.out
	}

	# generate VEP annotations
	call annotation.run_vep {
		input:
			ref = ref_ver,
			vcf = txt_to_vcf.out,
			vcf_idx = txt_to_vcf.idx,
			cache_dir = cache_dir
	}

	# Parse and append VEP columns to original vcf file
	call annotation.add_vep_cols {
		input:
			original_variants = gather_shards.out,
			vep_vcf = run_vep.vep_out
	}

	# flag GATK RankSum values (similar to samtools PV4)
	call annotation.flag_PV4 {
		input:
			infile = add_vep_cols.out
	}

	# flag SB (strand bias) 
	call annotation.flag_SB {
		input:
			infile = flag_PV4.out
	}

	# flag FDR (FDR-based min altdp) 
	call annotation.flag_FDR {
		input:
			infile = flag_SB.out
	}

	# flag RR (repeat region) 
	call annotation.flag_RR {
		input:
			infile = flag_FDR.out,
			map = rr_map,
			seg = rr_seg,
			lcr = rr_lcr
	}

	# flag VC (variant cluster) 
	call annotation.flag_VC {
		input:
			infile = flag_RR.out
	}

	#########################################################
	## Steps from Apply Filters
	#########################################################
	# estimate cohort AF from raw de novos file
	call annotation.estimate_cohort_AF {
		input:
			infile = flag_VC.out,
			cohort_size = cohort_size
	}

	# run CAF (cohort allele frequency)
	call annotation.flag_CAF {
		input:
			infile = flag_VC.out,
			caf_file = estimate_cohort_AF.out
	}

	#run outlier filter
	call annotation.flag_outlier {
		input:
			infile = flag_CAF.out,
			cohort_size = cohort_size,
			exp = expected_dnsnvs,
			cutoff = case_cutoff
	}

	#########################################################
	## parse filter flags, summarize filtering, output variants passing all filters
	#########################################################
	#run update_filter_column script to combine filter flags into single column
	call annotation.update_filt_col {
		input:
			infile = flag_outlier.out
	}

	#run script to summarize counts of variants flagged by each filter
	call annotation.summarize_counts {
		input:
			infile = update_filt_col.outfile
	}

	#run script to write out variants passing all filters, to be used as input to EM-mosaic
	call annotation.print_pass_vars {
		input:
			infile = update_filt_col.outfile
	}

	#########################################################
	## Run EM-mosaic
	#########################################################
	#run EM_mosaic
	call call_mosaic.detect_mosaic {
		input:
			infile = update_filt_col.outfile,
			postcut = postcut,
			outprefix = output_prefix
	}

	output {
		File denovos = detect_mosaic.denovos
		File mosaics = detect_mosaic.mosaics
		Array[File] output_plots = detect_mosaic.output_plots
	}


}