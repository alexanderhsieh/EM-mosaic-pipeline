## Usage: python flag_vclust.py <variants> <distance in bp> <output filename>
import sys

d = sys.argv[2] # distance to define a variant cluster

outf = open(sys.argv[3], 'w')

idx = {}

## add all variants to dict in form samples[id][chr] : [all variant positions]
allvars = []
samples = {} # id : [{chr1: [], chr2: [], etc}]
fullvars = {} # id|chr|pos : '\t'.join(tmp)
head = ''
with open(sys.argv[1], 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			head = '\t'.join(tmp) + '\t' + 'cluster_%s'%(d)
		else:
			allvars.append('\t'.join(tmp))
			id = tmp[0]
			if not id in samples:
				samples[id] = {}
			chr = tmp[idx['chr']]
			if not chr in samples[id]:
				samples[id][chr] = []
			pos = tmp[idx['pos']]
			samples[id][chr].append(pos)


#print head
outf.write(head + '\n')
for v in allvars:
	tmp = v.split('\t')
	if not tmp[0] == "id":
		id = tmp[idx['id']]
		chr = tmp[idx['chr']]
		pos = tmp[idx['pos']]
		# for each variant, flag as complex is it is located within <d> bp of another variant for same sample and chr
		complexflag = False
		if id in samples:
			if chr in samples[id]:
				for p in samples[id][chr]:
					dist = abs(int(pos) - int(p))
					if dist > 0 and dist <= int(d):
						complexflag = True
		#print '\t'.join(tmp) + '\t' + str(complexflag)
		out = '\t'.join(tmp) + '\t' + str(complexflag)
		outf.write(out + '\n')

f.close()

outf.close()