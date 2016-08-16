#!/usr/bin/python

import os
import sys
import random

def data_crossver(cv_path,dict_path,eventlist,cvlist,nfold):
	n = 0
	eventset = []
	nfold = int(nfold)
	names = locals()

        f = open(os.path.join(cv_path,cvlist),"w")

	for x in xrange(nfold):
		names['cvset%d' % x] = []
		names['cv_obj%d' % x] = open(os.path.join(cv_path,('cv%d' % x)),'w')
                print >>f, "cv%d" % x

        f.close()

	try:
		eventlist_obj = open(os.path.join(dict_path,eventlist),'r')

	except IOError as err:
		exit("You need run data_dict.py first")

	finally:
		for event in eventlist_obj.readlines():
			event = event.strip()
			wav_obj = open(os.path.join(dict_path,event),'r')
			for wav in wav_obj.readlines():
#				randomx = random.randint(0,4)
				wav = wav.strip()
				eventset.append([wav])
				n += 1
			random.shuffle(eventset)
			numstp = n/nfold
			for x in xrange(nfold):
				names['cvset%d' % x].append([event])
				for y in xrange(x*numstp,(x+1)*numstp):
					names['cvset%d' % x].append(eventset[y])
			print n
			eventset = []
			n = 0
		for x in xrange(nfold):
			for k in  names['cvset%d' % x]:
				print >> names['cv_obj%d' % x], k
			names['cv_obj%d' % x].close()

	print "****************************************"
	print "* Preparing Crossdata Completed!       *"
	print "****************************************"

if __name__ ==  '__main__':
	if len(sys.argv) != 6:
		exit('Usage: data_crossver.py <cv_path> <dict_path> <eventlist> <cvlist> <nfold>')
	data_crossver(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
