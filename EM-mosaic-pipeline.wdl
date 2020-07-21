version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines the workflow for using EM-mosaic to call mosaic SNVs from 
## single-sample gVCFs belonging to trios.  
## Note: currently limited to human whole-exome-sequencing data
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
		File parse_table_script
		File ref_fasta
		File ref_fasta_index
		File dn_script 
		Float pb_min_vaf
		Int par_min_dp
		Int par_max_alt
		String output_prefix
		String output_suffix

		## annotation
		String ref_ver
		File cache_dir
		File convert_script
		File parser_script
		File script_pv4 
		File script_sb 
		File script_fdr 
		File script_rr_parse 
		File rr_map 
		File rr_seg 
		File rr_lcr 
		File script_vc  

		## apply filters
		String output_prefix
		String CAF_outprefix
		File estimation_script
		Int cohort_size
		File script_CAF
		File outlier_script
		Int expected_dnsnvs 
		Int case_cutoff 
		File script_update_filter_col
		File script_filtct 
		File script_printpass 

		## filtered denovo to mosaic
		File em_mosaic_script
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
			script = parse_table_script,
			terra_table = sample_table,
			sample_map = sample_map,
			ped = ped
	}

	## Read in table and coerce as Array[File] corresponding to trios
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
				script = dn_script,
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
			variants = gather_shards.out,
			script = convert_script
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
			vep_vcf = run_vep.vep_out,
			script = parser_script
	}

	# flag GATK RankSum values (similar to samtools PV4)
	call annotation.flag_PV4 {
		input:
			infile = add_vep_cols.out,
			script = script_pv4
	}

	# flag SB (strand bias) 
	call annotation.flag_SB {
		input:
			infile = flag_PV4.out,
			script = script_sb
	}

	# flag FDR (FDR-based min altdp) 
	call annotation.flag_FDR {
		input:
			infile = flag_SB.out,
			script = script_fdr
	}

	# flag RR (repeat region) 
	call annotation.flag_RR {
		input:
			infile = flag_FDR.out,
			script_parse = script_rr_parse,
			map = rr_map,
			seg = rr_seg,
			lcr = rr_lcr
	}

	# flag VC (variant cluster) 
	call annotation.flag_VC {
		input:
			infile = flag_RR.out,
			script = script_vc
	}

	#########################################################
	## Steps from Apply Filters
	#########################################################
	# estimate cohort AF from raw de novos file
	call annotation.estimate_cohort_AF {
		input:
			script = estimation_script,
			infile = flag_VC.out,
			cohort_size = cohort_size
	}

	# run CAF (cohort allele frequency)
	call annotation.flag_CAF {
		input:
			infile = flag_VC.out,
			script = script_CAF,
			caf_file = estimate_cohort_AF.out
	}

	#run outlier filter
	call annotation.flag_outlier {
		input:
			infile = filter_CAF.out,
			script = outlier_script,
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
			infile = filter_outlier.out,
			script = script_update_filter_col
	}

	#run script to summarize counts of variants flagged by each filter
	call annotation.summarize_counts {
		input:
			infile = update_filt_col.outfile,
			script = script_filtct
	}

	#run script to write out variants passing all filters, to be used as input to EM-mosaic
	call annotation.print_pass_vars {
		input:
			infile = update_filt_col.outfile,
			script = script_printpass
	}

	#########################################################
	## Run EM-mosaic
	#########################################################
	#run EM_mosaic
	call call_mosaic.detect_mosaic {
		input:
			infile = update_filt_col.outfile,
			script = em_mosaic_script,
			postcut = postcut,
			outprefix = output_prefix
	}



	output {
		File denovos = detect_mosaic.denovos
		File mosaics = detect_mosaic.mosaics
		Array[File] output_plots = detect_mosaic.output_plots
	}


}