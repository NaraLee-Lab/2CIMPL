import os
import sys
import time
import re
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import math
import numpy as np
import scipy
import jpk_util
from matplotlib import figure
from PIL import Image
from PIL import ImageChops
from subprocess import call


wsn_chroms = jpk_util.getDir('../genome/WSN/WSN.chrom.sizes')

# Path used for annotation of the upper triangle plots
SEGMENT_AXES = ('/Users/jpk90/Desktop/scripts/python/img/segment_axes')


def getChromBins(chrom_len, bin_len):
	if bin_len > chrom_len:
		return [0]
	num_windows = int(chrom_len / bin_len)
	leftover_chrom = int(chrom_len % num_windows / 2)
	bins = [(leftover_chrom) + bin_len * i for i in range(1, num_windows)]
	bins = [0] + bins
	return bins


def getBins(bins, start, end, total_len):
	if start < 0:
		start = 0
		#print('Coordinates must be >= 0; got ' + str(start) + ', ' + str(end))
		#return None
	
	if end <= start:
		print('Error: end coordinate should never be before or equal to ' +
			  'the start coordinate')
		return None
	
	start_bin, end_bin = -1, -1
	start_diff, end_diff = 0, 0
	start_bin_size, end_bin_size = 0, 0
	
	for i, st in enumerate(bins):
		#print(str(i) + ': ' + str(st))
		if end_bin >= 0:
			break
		if start_bin < 0 and start < st:
			start_bin = i - 1
			start_bin_size = bins[i] - bins[i-1]
			start_diff = st - start
			#print('start_bin = ' + str(start_bin))
		if start_bin >= 0 and end < st:
			end_bin = i - 1
			end_bin_size = bins[i] - bins[i-1]
			end_diff = end - bins[i-1]
			#print('end_bin = ' + str(end_bin))
	
	if start_bin < 0:
		start_bin = len(bins) - 1
		start_bin_size = total_len - bins[-1]
		start_diff = total_len - start
	if end_bin < 0:
		end_bin = len(bins) - 1
		end_bin_size = total_len - bins[-1]
		end_diff = end - bins[-1]
			
	if start_bin != end_bin:
		
		#print('start_diff = ' + str(start_diff))
		#print('end_diff = ' + str(end_diff))
		
		#print('start_bin_size = ' + str(start_bin_size))
		#print('end_bin_size = ' + str(end_bin_size))
		
		read_25_pct = int((end - start) * 0.25) + 1
		if start_diff < read_25_pct and start_diff < (start_bin_size * 0.25):
			start_bin += 1
		if end_diff < read_25_pct and end_diff < (end_bin_size * 0.25):
			end_bin -= 1
		
		# check for whether or not the read should be in each bin
		
		# always round down + 1, so min 25% length = 1
		#read_25_pct = int((end - start) * 0.25) + 1
		#first_bin_pct = bins[start_bin + 1] - bins[start_bin]
		#if start_diff
	return (start_bin, end_bin)


def getBinnedLinkMatrix(link_file, bin_dict, lens_dict,
	include_nojuncs=True, just_inter=False, only_link_ends=False):
	'''
	include_nojuncs = if no junction is in the bin, still return it with a value of 0
	only_link_ends = junctions will only be considered at the end bins of each read
	'''
	junctions = {}
	if include_nojuncs:
		for chr1 in bin_dict:
			for chr2 in bin_dict:
				for i, b1 in enumerate(bin_dict[chr1]):
					for j, b2 in enumerate(bin_dict[chr2]):
						# only include if they're inter-segmental, unless the
						# "just_inter" flag is set to False
						if not just_inter or chr1 != chr2: 
							junc_name = chr1 + ':' + str(i) + '-' + chr2 + ':' + str(j)
							junc_name_rev = chr2 + ':' + str(j) + '-' + chr1 + ':' + str(i)

							#### New ####
							if junc_name_rev not in junctions:
								junctions[junc_name] = 0
							#############

							''' Old: may need to reuse
							if junc_name not in junctions:
								junctions[junc_name] = 0
							if junc_name_rev not in junctions:
								junctions[junc_name] = 0
							'''
	for i, line in enumerate(open(link_file, 'r')):
		#print(line.strip())
		chr1, st1, end1, chr2, st2, end2 = line.strip().split()[0:6]
		if chr1 not in bin_dict or chr2 not in bin_dict:
			print('Error on line ' + str(i) +
				  ': chromosome not found')
			return
		
		if just_inter and chr1 == chr2: 
			continue
		
		st1 = int(st1)
		st2 = int(st2)
		end1 = int(end1)
		end2 = int(end2)
		
		bins1 = getBins(bin_dict[chr1], st1, end1, lens_dict[chr1])
		bins2 = getBins(bin_dict[chr2], st2, end2, lens_dict[chr2])
		
		if bins1 == None or bins2 == None:
			continue

		if not only_link_ends:
			bins1 = range(bins1[0], bins1[1]+1)
			bins2 = range(bins2[0], bins2[1]+1)
		
		names = set()
		for b1 in bins1:
			for b2 in bins2:
				junc_name = chr1 + ':' + str(b1) + '-' + chr2 + ':' + str(b2)
				junc_name_rev = chr2 + ':' + str(b2) + '-' + chr1 + ':' + str(b1)
				names.add((junc_name, junc_name_rev))
				
		for k in names:

			#### New ####
			if k[0] in junctions:
				junctions[k[0]] += 1
			elif k[1] in junctions:
				junctions[k[1]] += 1
			else:
				junctions[k[0]] = 1
			#############

			''' Old
			junctions[k[0]] += 1
			if k[0] != k[1]:
				junctions[k[1]] += 1
			'''
	return junctions


def drawContactPlots(link_matrix, chrom_file=wsn_chroms,
	dir='.', prefix='', use_log_scale=False,
	cross_segment_intensities=False):

	max_val = max(link_matrix.values())
	if use_log_scale:
		max_val = math.log(max_val + 1)

	seg_lens = jpk_util.getChromSizes(chrom_file)
	for seg in seg_lens:
		pattern = re.compile(r'' + seg + ':[0-9]+-' + seg + ':[0-9]+')
		s = [(x,link_matrix[x]) for x in link_matrix.keys() if re.match(pattern, x)]

		# With new change, updated to be len * 2
		n = int(math.sqrt(len(s)*2))


		matrix = np.zeros(shape=(n,n))
		for i, info in enumerate(s):
			idx1 = info[0].index(':')
			idx2 = info[0].index('-')
			idx3 = idx2 + info[0][idx2:].index(':')

			c1 = int(info[0][idx1+1:idx2])
			c2 = int(info[0][idx3+1:])
			val = int(info[1])

			if val != 0:
				if use_log_scale:
					# add 1 so values of "1" don't register as 0
					val = math.log(1 + val)
				matrix[c1][c2] = val
				matrix[c2][c1] = val

		mask = np.tri(matrix.shape[0], k=-1)
		matrix = np.ma.array(matrix, mask=mask)

		fig, ax = plt.subplots()
		fig.set_size_inches(4, 4)
		ax.axis('off')

		if cross_segment_intensities:
			# For comparing intrasegments junctions from segment to segment,
			# rather than scaling the intensity to junctions within a segment
			ax.imshow(matrix, cmap='Blues', interpolation='nearest', vmax=max_val)
		else:
			ax.imshow(matrix, cmap='Blues', interpolation='nearest')

		image_name = dir + '/' + prefix + '_' + seg + '.png'
		if use_log_scale:
			image_name = dir + '/' + prefix + '_' + seg + '_logVals.png'
		if cross_segment_intensities:
			image_name = image_name.replace('.png', '_cross_intensities.png')

		plt.savefig(image_name, transparent=True, dpi=282)
		im = Image.open(image_name)
		im = im.rotate(45, expand=True)

		size = 1230, 999
		seg_axis = Image.open(SEGMENT_AXES + '/' + seg + '.png')
		seg_axis.thumbnail(size, Image.ANTIALIAS)
		im.paste(seg_axis, box=(198, 792))
		im.save(image_name)


def makeBinnedFile(link_matrix, bin_len, dir='.', prefix=''):
    bins_with_vals = [(x, link_matrix[x]) for x in link_matrix if link_matrix[x] > 0]
    dir = jpk_util.getDir(dir)

    binned_links_file = (dir + '/' + prefix + '_binned_' +
                         str(bin_len) + '_links.txt')
    out_f = open(binned_links_file, 'w+')
    for bin, numLinks in bins_with_vals:
        print(bin + "\t" + str(numLinks), file=out_f)
    out_f.close()


def generateContactPlots(links_file, bin_length,
	chrom_file=wsn_chroms, dir='.', prefix='', make_binned_file=True,
	use_log_scale=False, cross_segment_intensities=False):

	plot_dir = jpk_util.getDir(dir + '/contact_plots')
	link_dir = jpk_util.getDir(dir + '/links')

	# dictionary of segments & the start coordinate of each of their buckets
	seg_lens = jpk_util.getChromSizes(chrom_file)
	bin_dict = {}
	for segment in seg_lens:
		bin_dict[segment] = getChromBins(seg_lens[segment], bin_length)

	# Create bin matrices
	link_matrix = getBinnedLinkMatrix(links_file, bin_dict, seg_lens)

	if make_binned_file:
		makeBinnedFile(link_matrix, bin_length, dir=link_dir, prefix=prefix)

	drawContactPlots(link_matrix, chrom_file, dir=plot_dir, prefix=prefix,
		use_log_scale=use_log_scale,
		cross_segment_intensities=cross_segment_intensities)

	return link_matrix

