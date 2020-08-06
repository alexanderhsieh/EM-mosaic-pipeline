version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines tasks for running EM-mosaic
##
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

#Scores variants in dnSNVs file
task detect_mosaic {
	input {
		File infile
		Int postcut 
		String outprefix
	}

	command {
		Rscript /opt/EM-mosaic.v1.R ~{infile} ~{outprefix} ~{postcut} 
	}

	runtime {
		docker: "alexanderhsieh/em-mosaic-base:6Aug2020"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File denovos = "${outprefix}.denovo.txt"
		File mosaics = "${outprefix}.candidates.txt"
		Array[File] output_plots = glob("*.pdf")
	}
}
