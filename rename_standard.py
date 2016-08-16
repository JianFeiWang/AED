#!/usr/bin/python

import os
import sys

def rename_standard(prefix,suffix,path):
	srcfiles = os.listdir(path)
	srcfiles.sort()
	index = 1
	for srcfile in srcfiles:
		print srcfile
		suf = os.path.splitext(srcfile)[1]
		suf.lower()
		if suf == suffix:
			destfile = prefix + "%04d"%(index) + suffix
			os.rename(os.path.join(path,srcfile),os.path.join(path,destfile))
			index += 1
	print "**************************************"
	print " Standardized naming completed!       *"
	print "**************************************"
if __name__=='__main__':
	if len(sys.argv) != 4:
		exit("Usage:rename_standard.py <prefix> <suffix> <path>")
	rename_standard(sys.argv[1],sys.argv[2],sys.argv[3])
