version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines tasks for calling de novo SNVs from a table of single-sample gVCF paths
##
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

## Note: for ped, requires standard 6-column format at minimum:
## 		family ID, sample ID, father ID, mother ID, sex, affected status
## Note: for terra table, requires 2 columns: entity:sample_id, path to output_vcf_index
## Note: for sample map, requires 2 columns: sample_id, path to output_vcf
task read_terra_table {
	input {

		File terra_table
		File sample_map
		File ped
	}

	command {
		python /opt/parse_sample_table.py -i ~{terra_table} -m ~{sample_map} -p ~{ped} > "sample_attributes.tsv"
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:6Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "sample_attributes.tsv"
	}
}

## Reads in sample_attributes.no_header.tsv to enable coercion from 
## google bucket_path (String) to corresponding gvcf file (File)
## Note: requires specific 7-column format:
##		 sample_id, sample_gvcf, father_gvcf, mother_gvcf, 
##		 sample_gvcf_index, father_gvcf_index, mother_gvcf_index, 
task read_table {
	input{
		File table
	}

	command { 
		echo "reading table" 
	}

	runtime {
		docker: "ubuntu:latest"
		preemptible: 3
		maxRetries: 3
	}

	output{
		Array[Array[String]] out = read_tsv(table)
	}

}


## Runs bcftools merge on the trio gvcf array supplied in read_table()
## Parses read group ids from merged trio gvcf to be used in de novo calling
## Note: requires array of gvcf indices to be present as well
task merge_trio_gvcf {
	input{
		String sample_id

		Array[File] trio_gvcf_array
		Array[File] trio_gvcf_index_array

		File ref_fasta
		File ref_fasta_index

		Int disk_size = 100 
	}

	String outfname = "~{sample_id}.TRIO.g.vcf.gz"

	command {

		bcftools merge -g ~{ref_fasta} -l ~{write_lines(trio_gvcf_array)} -o ~{outfname} -O z -m all

		tabix -p vcf ~{outfname}

		bcftools query -l ~{outfname} > "ids.txt"


	}

	runtime {
		docker: "gatksv/sv-base-mini:cbb1fc"
		memory: "12G"
	    disks: "local-disk " + disk_size + " HDD"
	    preemptible: 3
		maxRetries: 3
	}

	output {
		File out_gvcf = "~{outfname}"
		File out_gvcf_index = "~{outfname}.tbi"
		Array[String] rg_ids = read_lines("ids.txt")
	}

}

## Calls denovos from bcftools merged gvcf containing proband, father, mother
## Note: readgroup IDs in the GVCF may differ from that in the PED or Sample Map
##		 here, we parse from the GVCF filename upstream when creating sample_attributes.tsv 
task call_denovos {
	input {
		String sample_id

		Array[String] trio_readgroup_ids

		File gvcf

		Float pb_min_vaf
		Int par_max_alt
		Int par_min_dp

		String output_suffix
	}

	String pb_id = trio_readgroup_ids[0]
	String fa_id = trio_readgroup_ids[1]
	String mo_id = trio_readgroup_ids[2]

	String output_file = "~{sample_id}~{output_suffix}"

	command {

		echo "~{pb_id} ~{fa_id} ~{mo_id}"

		python /opt/merged_gvcf_to_denovo.py -s ~{sample_id} -r ~{pb_id} -f ~{fa_id} -m ~{mo_id} -g ~{gvcf} -x ~{pb_min_vaf} -y ~{par_max_alt} -z ~{par_min_dp} -o ~{output_file}

		grep "^id" ~{output_file} > "header.txt"

		## fix id discordance
		#echo "converting ~{pb_id} to ~{sample_id}"

		#sed -i 's/~{pb_id}/~{sample_id}/g' ~{output_file}

	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:6Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{output_file}"
		File head = "header.txt"
	}
}

#Gathers shards of de novo calls into a single callset
task gather_shards {
	input {
		Array[File] shards 
		Array[File] headers
		String output_prefix
		String output_suffix
	}

	File header = headers[0]
	String output_file = "~{output_prefix}~{output_suffix}"

	command {

		while read file; do
		cat $file | grep -v "^id" >> "tmp.cat.txt"
		done < ~{write_lines(shards)};


		## CHECK THAT DE NOVO CALLSETS ARE NOT EMPTY
		N_LINES=`wc -l "tmp.cat.txt"`
		echo "tmp cat: $N_LINES lines"

		if [[ $(wc -l < "tmp.cat.txt") -le 1 ]]; then
			echo "EMPTY OUTPUT DE NOVO CALLSET"; exit $ERRCODE; 
		else
			(cat "~{header}" "tmp.cat.txt") > "~{output_file}"
		fi

		



	}

	runtime {
		docker: "ubuntu:latest"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{output_file}"
	}
}
