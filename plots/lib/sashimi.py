import os
import sys
import jpk_util

def generateSashimi(links_file, dir='.', prefix=''):
	'''
	links_file should be the output from makeLinkFile
	dir is the output directory
	prefix is the prefix for the sashimi bed file 
	'''
	dir = jpk_util.getDir(dir)

	sashimi_bed = dir + '/' + prefix + '_sashimi.bed'
		
	f1 = open(links_file, 'r')
	f2 = open(sashimi_bed, 'w+')

	# This will allow IGV to interpret the file as junctions
	f2.write('track name=%s_junctions or graphType=junctions\n' % prefix)

	for line in f1:
		ln = line.strip().split()
		if (len(ln) < 6):
			continue
			
		seg1 = ln[0]
		start1 = int(ln[1])
		end1 = int(ln[2])
		seg2 = ln[3]
		start2 = int(ln[4])
		end2 = int(ln[5])
		intra_segmental = seg1 == seg2
		more_than_3_nt = abs(start1 - start2) > 3
		non_overlapping = ((start1 < start2 or start1 > end2) and
			(start2 < start1 or start2 > end1))
		
		if intra_segmental and more_than_3_nt and non_overlapping:
			s1 = min(start1, start2)
			s2 = max(start1, start2)
			e1 = min(end1, end2)
			e2 = max(end1, end2)
			f2.write(('{}\t{}\t{}\t.\t1\t-\t{}\t{}\t150,50' +
				',50,30\t2\t{},{}\t0,{}\n').format(seg1,s1,e2,s1,e2,e1-s1,e2-s2,s2-s1))
	f1.close()
	f2.close()
