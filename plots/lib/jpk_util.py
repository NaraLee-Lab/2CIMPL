import os
import sys

# Colors for the circos interactions 
circos_color = {
	'PB2':'color=chr4_a3',
	'PB1':'color=chr6_a3',
	'PA':'color=chr8_a3',
	'HA':'color=chr17_a3',
	'NA':'color=chr12_a3',
	'NP':'color=chr14_a3',
	'M':'color=chr16_a3',
	'NS':'color=chr19_a3'
}

# Creates & gets the absolute path of the directory
def getDir(dir):
	dir = os.path.abspath(dir)
	if not os.path.exists(dir):
		try:
			os.mkdir(dir)
		except OSError:
			print('Directory was not found and could not be created: "' + dir + '"')
			exit(1)
	return dir


# To get a size dictionary for chromosomes
def getChromSizes(chrom_size_file):
	chrom_dict = {}
	for line in open(chrom_size_file, 'r'):
		chr, size = line.strip().split()
		chrom_dict[chr] = int(size)
	return chrom_dict


# Creates link files from Aligater Junction Files
def makeLinkFile(junction_file, dir='.', prefix='', out_name='',
	colored=False):
	'''
	junction_file should be the output from Aligater's "detect" function
	dir is the output directory
	prefix is the prefix for the links file 
	'''
	dir = getDir(dir)
	link_file = dir + '/' + prefix + '_links.txt'
	if (out_name != ''):
		link_file = dir + '/' + out_name

	out_f = open(link_file, 'w+')
	for line in open(junction_file, 'r'):
		ln = line.strip().split('\t')
		segs = ln[5].split(':')
		pos = ln[13].split(',')
		lens = ln[15].split(',')
		seqs = ln[10].split('_')
		read_id = ln[9]

		# usually there should only be two sections linked, but if there
		# are more, need to split them into separate records
		for i, seg1 in enumerate(segs):
			if len(segs) == i + 1:
				# reached the end
				break
			seg2 = segs[i + 1]
			flu1_len, flu2_len = lens[i: i + 2]
			flu1_seq, flu2_seq = seqs[i: i + 2]

			# For the position list, the first one is the end position of
			# the sequence while the rest are the start positions. I assume
			# this is because the ligation happened between the end position
			# of the first read and the next position of the second read,
			# though I'm not sure about after the second read since Aligater's
			# documentation says it should be going from that 2nd one to the
			# 3rd, I think
			if (i == 0):
				flu1_end, flu2_start = pos[i:i + 2]
				flu1_start = str(int(flu1_end) - int(flu1_len))
				flu2_end = str(int(flu2_start) + int(flu2_len))
			else:
				flu1_start, flu2_start = pos[i:i + 2]
				flu1_end = str(int(flu1_start) + int(flu1_len))
				flu2_end = str(int(flu2_start) + int(flu2_len))

			output = [seg1, flu1_start, flu1_end, seg2, flu2_start, flu2_end]
			if colored:
				output.append(circos_color[seg1])
			output_line = '\t'.join(output)

			print(output_line, file=out_f)
	out_f.close()
	return link_file


# Create a new link file from an existing one based on links that
# originate or end within a certain section
def makeLinkSubsetFile(link_file, chrom, chrom_size_file, start=-1,
	end=-1, dir='.', prefix='', out_name=''):
	dir = getDir(dir)
	chrom_lens = getChromSizes(chrom_size_file)

	if chrom not in chrom_lens:
		print('Chromosome %s was not found in chrom file' % chrom)
		exit(2)

	if start < 0:
		start = 0
	if end < 0:
		end = chrom_lens[chrom]
	if start > end:
		print('Invalid Coordinates; start (%s) should come before end (%s)'
			% (start, end))
		exit(2)

	subset_file = dir + '/' + prefix + '_subset.links'
	if (out_name != ''):
		subset_file = dir + '/' + out_name

	f_out = open(subset_file, 'w+')
	for line in open(link_file, 'r'):
		seg1, st1, end1, seg2, st2, end2 = line.strip().split()
		if ((seg1 == chrom and int(st1) > start and int(end1) < end) or
    		(seg2 == chrom and int(st2) > start and int(end2) < end)):
			print(line.strip(), file=f_out)
	f_out.close()
	return subset_file

