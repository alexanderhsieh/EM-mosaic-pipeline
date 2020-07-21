## Purpose: adds columns with GATK RankSum fields
## Usage: python filter_GATK_RankSum.py <VCF> <output filename>
## Output: writes to <output filename>
import sys
import scipy.stats as st


outf = open(sys.argv[2],'w')


with open(sys.argv[1],'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			head = '\t'.join(tmp) + '\t' + '\t'.join(['BaseQRankSum_z', 'MapQRankSum_z', 'ReadPosRankSum_z', 'BaseQRankSum_p', 'MapQRankSum_p', 'ReadPosRankSum_p'])
			outf.write(head + '\n')
		else:
			## parse INFO field and convert to dict
			info = tmp[idx['INFO']].split(';')
			keys = [i.split('=')[0] for i in info]
			vals = [i.split('=')[1] for i in info]
			infod = dict(zip(keys, vals))

			#print(infod)
			
			## parse key fields from INFO
			
			# initialize values
			mq_z, bq_z, rp_z = '.', '.', '.'

			if 'MQRankSum' in infod:
				mq_z = infod['MQRankSum']
			if 'BaseQRankSum' in infod:
				bq_z = infod['BaseQRankSum']
			if 'ReadPosRankSum' in infod:
				rp_z = infod['ReadPosRankSum']

			if not mq_z in ['.', 'NaN']:
				mq_p = st.norm.cdf(float(mq_z))
			if not bq_z in ['.', 'NaN']:
				bq_p = st.norm.cdf(float(bq_z))
			if not rp_z in ['.', 'NaN']:
				rp_p = st.norm.cdf(float(rp_z))

			#print([mq_z, bq_z, rp_z])
			#print([mq_p, bq_p, rp_p])


			
			out = '\t'.join(tmp) + '\t' + '\t'.join(map(str, [bq_z, mq_z, rp_z, bq_p, mq_p, rp_p]))
			outf.write(out + '\n')
outf.close()


