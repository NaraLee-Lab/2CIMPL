import os
import sys
import subprocess
import jpk_util
import contactPlots as cp
from subprocess import call


# Circos template configuration path
conf_template = jpk_util.getDir('../plots_lib/circos/template.conf')
influenza_karyotype = jpk_util.getDir('../plots_lib/circos/influenza_karyotype.txt')
wsn_chroms = jpk_util.getDir('../genome/WSN/WSN.chrom.sizes')

def makeIntersegFile(link_file, dir='.', prefix='', out_name=''):
	'''
	Creates a link file of just intersegmental ineractions
	'''
	dir = jpk_util.getDir(dir)
	interseg_file = dir + '/' + prefix + '_inter.links'
	if (out_name != ''):
		interseg_file = dir + '/' + out_name

	awk_cmd = 'awk \'$1!=$4\' ' + link_file + ' > ' + interseg_file
	os.system(awk_cmd)
	return interseg_file


def makeJunctionSiteFile(link_file, dir='.', prefix='', out_name='',
	num_nt_from_junc=10):
	'''
	Creates a link file of positions 10 nucleotides into the read from
	the junction sites. So if the link file says:
		PA 	10 	50 	NA 	7 	30
	Then the junction file will say:
		PA 	40 	50 	NA 	7 	17
	'''
	dir = jpk_util.getDir(dir)
	junc_file = dir + '/' + prefix + '_inter_junctions.links'
	if (out_name != ''):
		junc_file = dir + '/' + out_name

	f_out = open(junc_file, 'w+')
	for line in open(link_file, 'r'):

		fields = line.split()
		if len(fields) < 6:
			continue
		seg1, st1, end1, seg2, st2, end2 = fields[0:6]

		if seg1 != seg2: # only look at intersegmental junctions
			st1_upd = str(int(end1)-num_nt_from_junc)
			end2_upd = str(int(st2)+num_nt_from_junc)
			out_line = '\t'.join([seg1, st1_upd, end1, seg2, st2, end2_upd])
			if len(fields) > 6: # add circos color if it's there
				out_line += '\t' + fields[6]
			print(out_line, file=f_out)

	f_out.close()
	return junc_file


def makeBedgraphFromLinks(link_file, dir='.', prefix='', out_name='',
	chrom_size_file=wsn_chroms):
	'''
	Generate a bedgraph from the links. This can act as a histogram for
	the interactions that can be added to the border of the Circos plot
	'''
	dir = jpk_util.getDir(dir)
	bedgraph = dir + '/' + prefix + '_inter.bedgraph'
	if (out_name != ''):
		bedgraph = dir + '/' + out_name

	# This awk cmd will split the links into two separate lines, so:
	#	PA 	10 	100 	NA 	20 	30
	# Will become:
	# 	PA 	10 	100
	# 	NA 	20 	30
	awk_cmd = '\'{print $1\"\\t\"$2\"\\t\"$3\"\\n\"$4\"\\t\"$5\"\\t\"$6}\''

	# To leverage genomeCoverageBed, need to use system commands
	pipeline = ' '.join(['awk', awk_cmd, link_file, '|', 'sort',
		'-k1,1', '-k2,2n', '|', 'genomeCoverageBed', '-bga', '-i', 'stdin',
		'-g', chrom_size_file, '>', bedgraph])
	os.system(pipeline)
	return bedgraph


def makeInteractionFiles(link_file, dir='.', prefix='', chrom_size_file=wsn_chroms,
	fromJuncs=True):
	'''
	Generates all three files used for Circos
	'''
	dir = jpk_util.getDir(dir)
	try:
		if (fromJuncs):
			junc_file = makeJunctionSiteFile(link_file, dir=dir, prefix=prefix)
			bg = makeBedgraphFromLinks(junc_file, dir=dir, prefix=prefix,
				chrom_size_file=chrom_size_file)
			return (junc_file, bg)
		else:
			bg = makeBedgraphFromLinks(link_file, dir=dir, prefix=prefix,
				chrom_size_file=chrom_size_file)
			return (None, bg)
	except Exception as e:
		raise e
		exit(2)


def makeCircosConf(inter_links, prefix='', dir='.', colored=False,
	hitsclip_bg_path=None, hist_path=None, use_redux=True, out_name='',
	ribbon=False):
	'''
	colored defines whether the junctions will be colored or not. Nara
		has mentioned that it is probably more useful if it's not colored,
		since 
	hitsclip_bg_path is an optional parameter to include a bedgraph of
		HITS-CLIP data which will be seen around the plot
	'''
	dir = jpk_util.getDir(dir)

	if not colored:
		need_new_link_file = False
		for line in open(inter_links, 'r'):
			if len(line.strip().split()) > 6:
				need_new_link_file = True
			break
		if need_new_link_file:
			bname = '.'.join(os.path.basename(inter_links).split('.')[0:-1])
			new_link_file = inter_links.replace(bname, bname + '_nocolor')
			f_out = open(new_link_file, 'w+')
			for line in open(inter_links, 'r'):
				fields = line.strip().split()
				if len(fields) >= 6: 
					print('\t'.join(fields[0:6]), file=f_out)
			inter_links = new_link_file
			f_out.close()

	template = ''
	out_circos_name = dir + '/' + prefix + '.conf'
	template = conf_template
	
	f_in = open(template, 'r')
	f_out = open(out_circos_name, 'w+')
	for line in f_in:
		ln = line.strip().split()
		if len(ln) > 0 and ln[0] == 'links_file':
			f_out.write('links_file = ' + inter_links + '\n')
		elif len(ln) > 0 and ln[0] == 'hist_file' and hist_path:
			f_out.write('hist_file = ' + hist_path + '\n')
		elif len(ln) > 0 and ln[0] == 'hits_clip_file' and hitsclip_bg_path:
			f_out.write('hits_clip_file = ' + hitsclip_bg_path + '\n')
		elif len(ln) > 0 and ln[0] == 'ribbon':
			if ribbon:
				f_out.write('ribbon = yes\n')
			else:
				f_out.write('ribbon = no\n')
		elif len(ln) > 0 and ln[0] == 'karyotype':
			f_out.write('karyotype = ' + influenza_karyotype + '\n')
		else:
			f_out.write(line)
	f_in.close()
	f_out.close()
	return out_circos_name


def runCircosFromLinks(links_file, dir='.', prefix='', chrom_size_file=wsn_chroms,
	colored=False, hitsclip_bg_path=None, include_hist=True, ribbon=False, fromJuncs=True):

	dir = jpk_util.getDir(dir)
	circos_dir = dir + '/circos'
	conf_dir = dir + '/circos/conf'

	juncs, bg = makeInteractionFiles(links_file, dir=circos_dir,
		prefix=prefix, chrom_size_file=chrom_size_file, fromJuncs=fromJuncs)

	if not include_hist:
		bg = None

	conf_file = makeCircosConf(links_file, prefix=prefix, dir=conf_dir,
		hitsclip_bg_path=hitsclip_bg_path, hist_path=bg, ribbon=ribbon,
		colored=colored)

	pipeline = ' '.join(['circos', '-conf', conf_file, '-outputdir', circos_dir,
		'-outputfile', prefix + '_circos'])
	os.system(pipeline)


def generateThresholdLinks(inter_links_file, threshold=3, window=20,
	dir='.', prefix='', chroms_size_file=wsn_chroms, out_file=''):

	seg_lens = jpk_util.getChromSizes(chroms_size_file)
	bin_dict = {}
	for segment in seg_lens:
		bin_dict[segment] = cp.getChromBins(seg_lens[segment],
			window)
	new_dir = jpk_util.getDir(dir)
	if dir == '.':
		new_dir = (new_dir + '/binned_' + str(window) +
			'_junc_sites_t' + str(threshold))

	link_matrix = cp.getBinnedLinkMatrix(inter_links_file,
		bin_dict, seg_lens, just_inter=True, include_nojuncs=False)

	binned_thresh_file = (new_dir + '/' + prefix + '_binned_' +
		str(window) + '_ge' + str(threshold) + '_links.txt')
	if out_file != '':
		binned_thresh_file = out_file

	f_out = open(binned_thresh_file, 'w+')
	for link in link_matrix:
		if link_matrix[link] >= threshold:
			tup1, tup2 = link.split('-')

			seg1, bin1 = tup1.split(':')
			next_bin = int(bin1) + 1
			end1 = seg_lens[seg1]
			if len(bin_dict[seg1]) > next_bin:
				end1 = bin_dict[seg1][next_bin] - 1
			st1 = bin_dict[seg1][int(bin1)]

			seg2, bin2 = tup2.split(':')
			next_bin = int(bin2) + 1
			end2 = seg_lens[seg2]
			if len(bin_dict[seg2]) > next_bin:
				end2 = bin_dict[seg2][next_bin] - 1
			st2 = bin_dict[seg2][int(bin2)]

			out_line = '\t'.join([seg1, str(st1), str(end1),
				seg2, str(st2), str(end2)])
			print(out_line, file=f_out)
	f_out.close()
	return binned_thresh_file

