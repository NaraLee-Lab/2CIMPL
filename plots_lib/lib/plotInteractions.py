import os
import sys
import getopt
import textwrap
import time

import circos
import contactPlots
import sashimi
import jpk_util

'''
Produces Circos, Sashimi, and Upper-Triangle Plots based on Aligater's
output file. Specifically tailored for detecting RNA-RNA interactions
from influenza using hiCLIP or LIGR.
'''

DEFAULT_BIN_SIZE = 50

wsn_chroms = ('../../genome/WSN/WSN.chrom.sizes')

verbose = False


def displayUsage():
	print(textwrap.dedent("""usage=$ plotInteractions.py [-h] -i
	<juncs.txt> [options]

	where:

		-h  show this text

	required:
		-i  output junction file from Aligater

	optional:
		-C 	Generate Circos Plots
		-S 	Generate Sashimi Plot
		-U 	Generate Contact Plots

		-v  verbose (will print messages to stdout)
		-p 	prefix for the outbound files (defaults to basename of 
			the Aligater junction file)
		-o  output directory (defaults to current directory)
		-s  chromosome size file (default is WSN)

	circos parameters:
		-c 	colored output for Circos (default is it won't be colored)
		-H 	do not include histogram of interactions (default is True)
		-g 	use this HITS-CLIP bedgraph as well (default is None)

	contact plot parameters:
		-b 	bin size (default is 50)
		-l 	use a logarithmic scale (default is False)
		-x 	use the maximum number of interactions between ALL segments
			when defining the upper intensity limit, rather than just
			the max interactions within that segment. Compatible with -l.
			(default is False)
	"""))


def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hCSUxclvi:o:p:s:Hg:b:')
	except getopt.GetoptError:
		displayUsage()
		sys.exit(2)

	doCircos = False
	doSashimi = False
	doContact = False

	junction_file = ''
	prefix = ''
	output = '.'
	chrom_size_file = wsn_chroms
	hits_clip_file = None
	bin_size = DEFAULT_BIN_SIZE
	use_log_scale = False

	colored = False
	include_hist = True
	cross_segment_intensities = False

	for opt, arg in opts:
		if opt == '-h':
			displayUsage()
			sys.exit()
		elif opt == "-i":
			junction_file = arg
		elif opt == "-o":
			output = arg
		elif opt == "-p":
			prefix = arg
		elif opt == "-v":
			verbose = True
		elif opt == "-c":
			colored = True
		elif opt == "-s":
			chrom_size_file = arg
		elif opt == '-H':
			include_hist = False
		elif opt == '-g':
			hits_clip_file = arg
		elif opt == '-b':
			bin_size = int(arg)
		elif opt == '-l':
			use_log_scale = True
		elif opt == '-x':
			cross_segment_intensities = True
		elif opt == '-C':
			doCircos = True
		elif opt == '-S':
			doSashimi = True
		elif opt == '-U':
			doContact = True

	chrom_size_file = jpk_util.getDir(chrom_size_file)

	if junction_file == '':
		displayUsage()
		sys.exit(2)

	if prefix == '':
		prefix = os.path.basename(junction_file)
		prefix = prefix[0:prefix.index('.')]

	links_dir = output + '/links'
	link_file = jpk_util.makeLinkFile(junction_file, dir=links_dir, prefix=prefix)

	if doCircos:
		inter_links = circos.makeIntersegFile(link_file, dir=output + '/circos',
			prefix=prefix)
		circos.runCircosFromLinks(inter_links, dir=output, prefix=prefix,
			chrom_size_file=chrom_size_file, colored=colored,
			hitsclip_bg_path=hits_clip_file, include_hist=include_hist)

	if doSashimi:
		sashimi.generateSashimi(link_file, dir=output + '/sashimi', prefix=prefix)

	if doContact:
		contactPlots.generateContactPlots(link_file, bin_size,
			chrom_file=chrom_size_file, dir=output, prefix=prefix,
			use_log_scale=use_log_scale,
			cross_segment_intensities=cross_segment_intensities)


if __name__ == '__main__':
	main()
