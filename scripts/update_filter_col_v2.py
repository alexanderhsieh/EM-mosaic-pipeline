## Purpose: parses all filter-related columns and creates an aggregated 'filter' column
## 			requires preprocessing PV4, SB, FDR, RR, VC
##			adds MAF during this step
## 			(optional/TBD) cohort AF (CAF) parsing
## NOTE: Assumes 'filter' column already exists (from reformat_vcf.py step upstream)
'''
Usage: python update_filter_col.py [options]

Output: writes to output_file

Options:
  -h, --help            show this help message and exit
  -i input_file, --input=input_file
                        input tab-separated variants file

  -o output_file, --output=output_file
'''
import sys
from optparse import OptionParser

########################################
## Handle arguments
########################################

parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-o', '--output', dest='output_file',help='output filename')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	sys.exit()

input_file = options.input_file
output_file = options.output_file

########################################
## create dictionary of cutoff values
########################################
#cutoff = {'PV4_bq': '1e-3', 'PV4_mq': '1e-6', 'PV4_rp': '1e-3', 'MAF': '1e-4', 'CAF': '0.01'}
cutoff = {'PV4_bq': '0.05', 'PV4_mq': '0.05', 'PV4_rp': '0.05', 'MAF': '1e-4', 'CAF': '0.003'}
  
#print('')
#print(cutoff)
#print('')
########################################
## Update variants file
########################################

filtf = open(output_file, 'w')
with open(input_file, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			if 'filter' in idx:
				filtf.write('\t'.join(tmp) + '\n')
			else:
				filtf.write('\t'.join(tmp) + '\t' + 'filter' + '\n')

			#print(idx)

		else:
			# initialize final filter
			outfilt = []

			# fetch relevant values
			pv4_bq, pv4_mq, pv4_rp = tmp[idx['BaseQRankSum_p']], tmp[idx['MapQRankSum_p']], tmp[idx['ReadPosRankSum_p']]
			sbflag = tmp[idx['strand_bias_flag']] # 0/1
			fdrmin = tmp[idx['fdr.min.alt']] 
			bed_map, bed_seg, bed_lcr = tmp[idx['map_score']], tmp[idx['segdup_score']], tmp[idx['LCR_flag']]
			vclust = str(tmp[idx['cluster_10']]) # TRUE / FALSE
			popfreq = tmp[idx['MAX_AF']]
			biotype = tmp[idx['BIOTYPE']]
			#refgene_func = tmp[idx['Func.refGene']] # DEPRECATED - find alternative
			gene = tmp[idx['SYMBOL']]
			#rsid = tmp[idx['rs_dbSNP151']] # DEPRECATED
			rsid = tmp[idx['Existing_variation']] # DEPRECATED

			share_flag = tmp[idx['shared_pb_sib']]
			

			if tmp[idx['cohort_AF']] in ['.', 'NA']:
				cohort_af = 0.0
			else:
				cohort_af = float(tmp[idx['cohort_AF']])

			if 'outlier_flag' in idx:
				outlier_flag = tmp[idx['outlier_flag']]
			else:
				outlier_flag = 'FALSE'
			if outlier_flag == '.':
				outlier_flag = 'FALSE'


			## filter logic
			
			## PV4
			if float(pv4_bq) < float(cutoff['PV4_bq']) or float(pv4_mq) < float(cutoff['PV4_mq']) or float(pv4_rp) < float(cutoff['PV4_rp']):
				outfilt.append('PV4_FAIL')
			else:
				outfilt.append('PV4_PASS')

			## Strand Bias
			if sbflag == '1':
				outfilt.append('SB_FAIL')
			else:
				outfilt.append('SB_PASS')

			## FDR-based minimum Nalt
			if float(tmp[idx['altdp']]) < float(fdrmin):
				outfilt.append('FDR_FAIL')
			else:
				outfilt.append('FDR_PASS')

			## Repeat Regions (mappability, segdup, LCR)
			if bed_map in ['.', 'NA'] and bed_seg in ['.', 'NA'] and bed_lcr == 'no':
				outfilt.append('RR_PASS')
			else:
				outfilt.append('RR_FAIL')

			## Variant cluster
			if vclust in ['True', 'TRUE', 'true']:
				outfilt.append('VC_FAIL')
			else:
				outfilt.append('VC_PASS')

			## Population Frequency
			if popfreq in ['NA', '.']:
				popfreq = 0.0
			if float(popfreq) > float(cutoff['MAF']):
				outfilt.append('MAF_FAIL')
			else:
				outfilt.append('MAF_PASS')

			## Protein-coding filter
			#if biotype == 'protein_coding' and refgene_func in ['exonic', 'splicing']:
			if biotype == 'protein_coding': # edit 04/16
				outfilt.append('COD_PASS')
			else:
				outfilt.append('COD_FAIL')

			## MUC/HLA gene filter
			if gene.startswith('MUC') or gene.startswith('HLA'):
				outfilt.append('MUC-HLA_FAIL')
			else:
				outfilt.append('MUC-HLA_PASS')

			## dbSNP filter
			if rsid.startswith('rs'):
				outfilt.append('dbSNP_FAIL')
			else:
				outfilt.append('dbSNP_PASS')

			## cohortAF
			if cohort_af >= float(cutoff['CAF']):
				#print(('FAIL', cohort_af, cutoff['CAF']))
				outfilt.append('CAF_FAIL')
			else:
				#print(('PASS', cohort_af, cutoff['CAF']))
				outfilt.append('CAF_PASS')

			if outlier_flag == 'TRUE':
				outfilt.append('OUT_FAIL')
			else:
				outfilt.append('OUT_PASS')

			if share_flag == 'True':
				outfilt.append('SHARE_FAIL')
			else:
				outfilt.append('SHARE_PASS')


			## update filter column
			## in the future, may skip having the filter column and just append at the end in this step...
			if 'filter' in idx:
				tmp[idx['filter']] = '|'.join(outfilt)
				## write out updated file
				tmpout = '\t'.join(tmp) 
				filtf.write(tmpout + '\n')

				#print(tmpout)
			
			else:
				tmpout = '\t'.join(tmp) + '\t' + '|'.join(outfilt)
				filtf.write(tmpout + '\n')

				#print(tmpout)


filtf.close()


