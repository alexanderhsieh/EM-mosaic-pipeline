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

import 'https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_gvcf_to_denovo.wdl' as gvcf_to_denovo

###########################################################################
# WORKFLOW DEFINITION
###########################################################################
workflow EM_mosaic_pipeline {

	input {

		File sample_table
		File sample_map
		File ped


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
	## Call de novo SNVs from gVCFs in sample map
	##########################################################################
	






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