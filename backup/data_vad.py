#!/usr/bin/python

import sys
import os
import numpy as np
#import multiprocessing as mp
from package.python_speech_features import fbank

def read_data(dict_path,eventlist):

	eventlist_obj = open(os.path.join(dict_path,eventlist))
	mevents = 0
	event_dict = {}

	for event in eventlist_obj.readlines():
		nevent = 0
		event = event.strip()
		print '%s' % event

		data_sum = []
		wavdata_obj = open(os.path.join(dict_path,event))

		for data in wavdata_obj.readlines():
#			print nevent
			data = data.strip()
			data_sum.append(data)
			nevent += 1
#			os.system("pause")
#		raw_input()
#		print data_sum
		event_dict[event] = data_sum
		wavdata_obj.close()
		mevents += nevent

	eventlist_obj.close()

	print "The number of events %d" % mevents

	return mevents,event_dict

def medfilt(data):

	for x in xrange(len(data)-2):
		data[x] = np.dot(data[x:x+3],np.ones(3).T)/3

	return data[:-2]

def time_feature_extract(wdata_frame,win_ham):

	Rxx = np.zeros(128)
	h = np.zeros(1)
	D = np.zeros(128)
	c = np.zeros(1)

	for k in xrange(128):
 		Rxx[k] = sum(wdata_frame[0:128] * wdata_frame[k:128+k])/sum(win_ham[0:128]*win_ham[k:k+128])

# feature extracted:
# A) Time Domain Features
#	clearity
	h = max(Rxx)/(Rxx[0]-max(Rxx[16:]))
#	harmonicity
	D = 0.8 * np.sqrt(abs(2*(Rxx[0]*np.ones(128) - Rxx)))
	c = 1 - min(D[16:])/max(D[16:])
#	predicition gain
	p = 10
	a = np.zeros(p)
	f = np.zeros(p)
	b = np.zeros(p)
	a[0] = -Rxx[1]/Rxx[0]
	f[0] = 1/Rxx[0]
	b[p-1] = 1/Rxx[0]

	for x in xrange(2,p):
		Ef = sum(Rxx[x-1:0:-1]*f[:x-1])
		Eb = sum(Rxx[1:x]*b[1-x:])
		Ea = sum(Rxx[x-1:0:-1]*a[:x-1])
		f1 = 1/(Eb*Ef)*f[:x]-Ef/(1-Eb*Ef)*b[-x:]
		b[-x:] = 1/(Eb*Ef)*b[-x:]-Eb/(1-Eb*Ef)*f[:x]
		f[:x] = f1
		a[:x] = a[:x] + (Rxx[x]-Ea)*b[-x:]

	Gp = np.log(abs(Rxx[0]/Ea))

	return h,c,Gp

def frequency_feature_extract(wdata_frame,wdata_frame2):

# B) Frequency Domain Features
	Xdata_frame_dft = np.fft.fft(wdata_frame,n=2048)
#	Xdata_frame_dft2 = np.fft.fft(wdata_frame2,n=2048)
#	Periodictity
	Per = np.zeros(1500)

	for w in xrange(60,1560):
		Per[w-60] =  np.log(abs(Xdata_frame_dft[int(32 / 125 * w)])) + np.log(abs(Xdata_frame_dft[int(2 * 32 / 125 * w)]))

	Phps = max(Per)

#	spectral flux
	feat1,enery1 = fbank(wdata_frame, samplerate = 8000, winstep = 0.032, nfilt = 80, nfft = 2048, lowfreq = 60)
	feat2,enery2 = fbank(wdata_frame2, samplerate = 8000, winstep = 0.032, nfilt = 80, nfft = 2048, lowfreq = 60)

	SFp = sum(abs(feat2[0] - feat1[0]))

	return Phps,SFp

def vad(dict_path,feat_path,eventlist,vadfile):

	Vad_vector = {}
	win_ham = np.hamming(256)
	Count = 0

	mevents,event_dict = read_data(dict_path,eventlist)

## 32ms hamming 10ms skip
#  compute the autocorrelation

	for key in event_dict:
		Vad_vector[key] = []

	print 'loading...'

	f = open(os.path.join(feat_path,vadfile),"w")

	for key in event_dict:
#	for key in ["sli","sil"]:

		ldata = event_dict[key]

		for data in ldata:

			Count += 1
			adata = []

			for sdata in data.split():
				adata.append(int(sdata))

			npoint =  len(adata)
			nframe = (npoint - (32 * 8)) / (10 * 8) + 1
			if nframe <= 0:
				continue

			data_frame = np.zeros((nframe,256))
			vector = np.zeros((nframe-1,5))

			for n in xrange(nframe):

				data_frame[n] = adata[(n * 80) : (256 + n * 80) ]

			for n in xrange(nframe-1):

				wdata_frame = data_frame[n] * win_ham
				wdata_frame2 = data_frame[n+1] * win_ham
#       multiprocessing	****************************#
#				if multi:
#					pool = mp.Pool(processes = 2)
#					h,c,Gp = pool.apply(time_feature_extract,(wdata_frame,win_ham,))
#					Phps,SFp = pool.apply(frequency_feature_extract,(wdata_frame,wdata_frame2,))
#					pool.close()
#					pool.join()
				#name = mp.current_process().name
				#print name, 'Starting'

				#worker_1 = mp.Process(name = 'worker_1', target = time_feature_extract,args=(wdata_frame,win_ham))
				#worker_2 = mp.Process(name = 'worker_2', target = frequency_feature_extract,args=(wdata_frame,wdata_frame2))

				#worker_1.start()
				#worker_2.start()
				#print name,"Exiting"
# ********************************************************
#				else:
				h,c,Gp = time_feature_extract(wdata_frame,win_ham)
				Phps,SFp = frequency_feature_extract(wdata_frame,wdata_frame2)
				vector[n] = [h,c,Gp,Phps,SFp]

			sys.stdout.write('.')
			sys.stdout.flush()

			if (Count%13==0):
				print " %.1f %%" % (Count/float(mevents) * 100)
# C) Post-handling
#	 normalization
			vector = (vector - np.mean(vector,axis = 0) )/(np.std(vector,axis = 0)+10E-5)
#	 PCA
			cov_vector = np.cov(vector.T)
			try:
				u,s,v = np.linalg.svd(cov_vector)
			except:
#				print vector
#				print cov_vector
				continue
			M = v[0]
			svd_vector = np.dot(vector,M.T)
#smooth
			Vad_vector[key].append(medfilt(svd_vector))

	for key,value in Vad_vector.items():
		print >>f, "%s : %s" % (key,value)
	f.close()
	print "100%"

	print "******************************************"
	print "* soft-decision prepared!                *"
	print "******************************************"
#hard-decision:
#	 2-class GMM







if __name__=="__main__":
	if len(sys.argv) != 5:
		exit("Usage: vad.py <dict_path> <feature_path> <eventlist> <vadvector> <multiprocess>")
	vad(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
