#!/usr/bin/python3
## Purpose: call de novos from a trio gvcf
'''
Usage: gvcf_to_denovo.py -s <sample id> \
                         -g <gzipped trio gvcf>
                         -p <PED file>
                         -x <proband min vaf> \
                         -y <parent max altdp> \
                         -z <parent min dp> \
                         -o <output filename> 

## CAVEATS:
# -SNVs only, no indels or SVs
# -splits multiallelic sites into individual lines
# -ignores missing genotypes ('./.') and sites without AD or DP information

# Output format:
# chr, pos, ref, alt, refdp, altdp, dp, adfref, adfalt, adrref, adralt, <proband VCF information repeated>,
#                <father FORMAT field>, <father GT information>,
#                <mother FORMAT field>, <mother GT information>
'''
import sys
from optparse import OptionParser
import subprocess
import gzip
import io
import os



####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-s', '--sid', dest='sample_id',help='sample id')
parser.add_option('-r', '--rgid', dest='rg_id',help='sample readgroup id')
parser.add_option('-f', '--faid', dest='faid', help='father id')
parser.add_option('-m', '--moid', dest='moid', help='mother id')
parser.add_option('-g', '--gvcf', dest='gvcf', help='trio gvcf')
parser.add_option('-x', '--min_vaf', dest='pb_min_vaf',help='proband minimum variant allele frequency')
parser.add_option('-y', '--max_alt', dest='par_max_alt',help='parent maximum alternate allele read depth')
parser.add_option('-z', '--min_dp', dest='par_min_dp',help='parent minimum read depth')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_id == None or options.rg_id == None or options.gvcf == None or options.faid == None or options.moid == None  or options.pb_min_vaf == None or options.par_max_alt == None or options.par_min_dp == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


sample_id = options.sample_id
rg_id = options.rg_id
faid = options.faid
moid = options.moid
gvcf = options.gvcf
pb_min_vaf = options.pb_min_vaf
par_max_alt = options.par_max_alt
par_min_dp = options.par_min_dp
output_file = options.output_file


####################################################################################################
## iterate over merged_trio gVCF and call denovos
####################################################################################################
#bufsize=1
#outf = open(output_file, 'w', buffering=bufsize)
outf = open(output_file, 'w')

cmd2 = 'zcat < %s | grep -v "#"| wc -l'%(gvcf)
#cmd2 = 'cat %s | grep -v "#"| wc -l'%(gvcf)
tot = int(subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip().split(' ')[0])
print('## TOTAL VARIANT LINES: %s'%(str(tot)))

print('')
print('## ITERATING OVER VARIANT LINES')
print('')

head = ['id', 'rg_id', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'dp', 'adfref', 'adfalt', 'adrref', 'adralt', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'S_GT', 'FA_GT', 'MO_GT']
#print '\t'.join(head)
outf.write('\t'.join(head) + '\n')

i = 0
dnct = 0


## iterate over gVCF
#with open(gvcf, 'r') as f:
#  for line in f:

with gzip.open(gvcf, 'rb') as f:
  with io.TextIOWrapper(f, encoding='utf-8') as decodef:  
    for line in decodef:


      tmp = line.strip().split('\t')

      ## handle vcf header information
      if line.startswith('##'):
        continue
      
      ## handle vcf header containing column information
      if line.startswith('#CHROM'):
        idx = {col:index for index, col in enumerate(tmp)}

      ## handle variant lines
      else:
        i += 1
        if i%100000 == 0:
          print('## %d/%d lines processed ... '%(i, tot))


        ## initialize values to avoid iteration bugs
        chr, pos, ref, alt = '', '', '',''
        region = ''
        info = ''
        fmt, gt = [], []
        gtd = {}


        ## ignore non-variant blocks 
        if not tmp[idx['INFO']].startswith('END='):
          
          ## get proband variant information
          chr, pos, ref = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']]

          # get region for tabixing parents
          region = chr + ':' + pos + '-' + pos



          ## parse alt allele
          #alt = tmp[idx['ALT']].strip(',<NON_REF>')
          alt = tmp[idx['ALT']].split(',')

          
          ## how to handle multiallelic sites? e.g. chr1    1646352 .       A       C,G,<NON_REF>
          ## iterate over all alternate alleles present in ALT
          #for a in alt.split(','):
          for a in alt:
            if not a == '<NON_REF>':
              if not a == '*': # ignore point deletions for now; messy when matching alleles with parents
                #if len(ref) == 1 and len(a) == 1: # ignore indels for now;
                if len(ref) == 1 and len(ref) == len(a):  

                  if len(tmp) == len(idx.keys()): # ignore lines with missing fields
                    
                    # Get index of current alternate allele
                    #pb_altidx = alt.split(',').index(a) + 1
                    pb_altidx = alt.index(a) + 1 # need to leave <NON_REF> in to correctly parse alt allele index


                    # save INFO field
                    info = tmp[idx['INFO']]

                    # create dictionary of FORMAT:GT mapping
                    fmt = tmp[idx['FORMAT']].split(':')
                    

                    #pb_gt = tmp[idx[sample_id]].split(':') ## ASSUMES THAT SAMPLE GENOTYPE INFORMATION IS IN THE LAST COLUMN; didn't use ID since column ID differs from sample id....
                    pb_gt = tmp[idx[rg_id]].split(':') ## ASSUMES THAT SAMPLE GENOTYPE INFORMATION IS IN THE LAST COLUMN; didn't use ID since column ID differs from sample id....
                    pb_gtd = dict(zip(fmt, pb_gt)) # e.g. {'GT': '0/1', 'AD': '5,7,0', 'GQ': '99', 'PL': '157,0,104,172,125,297', 'SB': '5,0,7,0', 'DP': '12'}

                    fa_gt = tmp[idx[faid]].split(':')
                    fa_gtd = dict(zip(fmt, fa_gt))

                    mo_gt = tmp[idx[moid]].split(':')
                    mo_gtd = dict(zip(fmt, mo_gt))

                    ## CHECK THAT ALL NECESSARY INFORMATION IS PRESENT
                    #if not (pb_gtd['GT'] == './.'): ## ignore sites with missing genotypes
                    #  if ('AD'in pb_gtd) and ('DP' in pb_gtd): # ignore sites with no AD or DP information
                    if not './.' in [pb_gtd['GT'], fa_gtd['GT'], mo_gtd['GT']]:
                      
                      if 'AD' in pb_gtd and 'DP' in pb_gtd and 'AD' in fa_gtd and 'DP' in fa_gtd and 'AD' in mo_gtd and 'DP' in mo_gtd:

                        if len(pb_gtd['AD'].split(',')) > 1: # ensure there is alternate allele read support
                          #if not pb_gtd['GT'] in ['0/0', '0|0']:
                            
                          #print(line)


                          ## parse strand-specific allelic depth information
                          adf = pb_gtd['F1R2']
                          adr = pb_gtd['F2R1']

                          adfref = adf.split(',')[0]
                          adrref = adr.split(',')[0]

                          adfalt = adf.split(',')[pb_altidx]
                          adralt = adr.split(',')[pb_altidx]

                          pb_refdp = pb_gtd['AD'].split(',')[0]
                          pb_altdp = pb_gtd['AD'].split(',')[pb_altidx]
                          if pb_altdp == '.':
                            pb_altdp = '0'
                          #pb_dp = pb_gtd['DP']
                          pb_dp = int(pb_refdp) + int(pb_altdp)

                          if int(pb_dp) > 0:
                            pb_vaf = float(pb_altdp)/float(pb_dp)
                          else:
                            pb_vaf = 0.0


                          ## parse parental information
                          fa_dp = fa_gtd['DP']
                          fa_refdp = fa_gtd['AD'].split(',')[0]

                          if fa_gtd['AD'] == '.': # handle case where AD is just .
                            fa_altdp = '0'
                          else:
                            fa_altdp = fa_gtd['AD'].split(',')[pb_altidx]
                            if fa_altdp == '.': # handle case where AD is 10,.,3
                              fa_altdp = '0'
                          
                          mo_dp = mo_gtd['DP']
                          mo_refdp = mo_gtd['AD'].split(',')[0]
                          #mo_altdp = mo_gtd['AD'].split(',')[pb_altidx]

                          if mo_gtd['AD'] == '.': # handle case where AD is just .
                            mo_altdp = '0'
                          else:
                            mo_altdp = mo_gtd['AD'].split(',')[pb_altidx]
                            if mo_altdp == '.': # handle case where AD is 10,.,3
                              mo_altdp = '0'
                          

                          ## APPLY DE NOVO CALLING CRITERIA
                          if not (int(pb_refdp) == 0): # ignore hom alt sites
                            if float(pb_vaf) >= float(pb_min_vaf):
                              if int(fa_altdp) <= int(par_max_alt) and int(mo_altdp) <= int(par_max_alt):
                                if int(fa_dp) >= int(par_min_dp) and int(mo_dp) >= int(par_min_dp):
                                  out = map(str, [sample_id, rg_id, chr.strip('chr'), pos, ref, a, pb_refdp, pb_altdp, pb_dp, adfref, adfalt, adrref, adralt])

                                  #print '\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1]
                                  #outf.write('\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1] + '\n')
                                  outstring = '\t'.join(out) + '\t' + '\t'.join(tmp[:idx['FORMAT']+1]) + '\t' + ':'.join(pb_gt) + '\t' +  ':'.join(fa_gt) + '\t' + ':'.join(mo_gt)
                                  #print(outstring)
                                  outf.write(outstring + '\n')
                                  # deal with empty output?
                                  outf.flush()
                                  os.fsync(outf)

                                  dnct += 1

                                  print('## %d de novo variants found ...'%(dnct))





outf.close()


