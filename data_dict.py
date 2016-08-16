#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import re
import wave
import numpy


def read_wav(wav_path,name):
	'''
		Usage:[SampleNumber, SampleRate, WaveDate] read_wav.py <wav_path> <name>
	'''
	datawav = wave.open(os.path.join(wav_path,name),'rb')
#	read paramters
#	return a turpler: track number,quantitative figures,sampling rate,frame length
	params = datawav.getparams()
	nchannels,sampwidth,framerate,nframes = params[:4]
	print "******************************************"
	print "* THE PARAMETERS OF THE AUDIO            *"
	print "* Channels : %-27d *" % nchannels
	print "* Quantity(byte) : %-21d *" % sampwidth
	print "* SampleRate : %-25d *" % framerate
	print "* SampleNumber : %-23d *" % nframes
	print "******************************************"
#	read a fixed length wav data
	str_data = datawav.readframes(nframes)
	datawav.close()
#	the wavform is converted into an array
	wave_data = numpy.fromstring(str_data,dtype = numpy.short)
#	wave_data.shape = -1,2
#	wave_data = wave_data.T
#	time = numpy.arange(0,nframes)*(1.0/framerate)
#	print time
	return nframes,framerate,wave_data

def read_textgrid(label_path,name):
	datalabel = open(os.path.join(label_path,name),'r')
	bPrint = 0
	minlist = []
	maxlist = []
	textlist = []
	for index in datalabel.readlines():
		if re.search('item \[1\]',index):
			bPrint = 1
		if bPrint:
			m = re.search(r'xmin',index)
			n = re.search(r'xmax',index)
			l = re.search(r'text',index)
			if m:
				v = re.findall(r'xmin = (\d+\.?\d*)',index)[0]
				xmin = int(float(v) * 8000)
			elif n:
				v = re.findall(r'xmax = (\d+\.?\d*)',index)[0]
				xmax = int(float(v) * 8000)
			elif l:
				v = re.findall(r'text = \"(.*)\"',index)[0]
				text = v.strip()
				if not text:
					text = 'oth'
				minlist.append(xmin)
				maxlist.append(xmax)
				textlist.append(text)
	return minlist,maxlist,textlist

def data_dict(path,wavlist,eventlist):
	WavNum = 0
	error = 0
	Wav_List = []
	Event_List = []
	try:
		wav_list_obj = open(os.path.join(path,wavlist),'r')
	except IOError as err:
		print "No find the list, we will make new one"
		error = 1
		info = os.path.join(path,'wav')
		listfile = os.listdir(info)
		try:
			wav_list_obj_b = open(os.path.join(path,wavlist),'w')
		except IOError as err:
			exit(err)
		for line in listfile:
			new_line = line
			print new_line
			new_line.lower()
			if new_line[-4:]=='.wav':
 				WavNum +=  1
				Wav_List.append(line[:-4])
		Wav_List.sort()
		wav_list_obj_b.write("The number of wavfile : "+str(WavNum)+'\n')
		for line in Wav_List:
			wav_list_obj_b.write(line+"\n")
		wav_list_obj_b.close()
		print "completed!"
	finally:
		if error:
			wav_list_obj = open(os.path.join(path,wavlist),'r')
		FL = wav_list_obj.readline()
		print FL
		WavNum = int(re.findall(r'(\d+)',FL)[0])
		wav_path = os.path.join(path,'wav')
		label_path = os.path.join(path,'label')
		dict_path = os.path.join(path,'dict')
		wav_dict = {}
		label_dict = {}
		v = 0
		for index in wav_list_obj.readlines():
			index = index.strip()
			wav = index + '.wav'
			label = index + '.TextGrid'
			try:
				v += 1
				print v
 				wav_pdata = read_wav(wav_path,wav)
 				wav_dict[index] = wav_pdata
				label_pdata = read_textgrid(label_path,label)
				label_dict[index] = label_pdata
			except TypeError as err:
 				wav_list_obj.close()
 				exit(err)

		wav_list_obj.seek(28)
		for index in wav_list_obj.readlines():
			index = index.strip()
			xmin = label_dict[index][0]
			xmax = label_dict[index][1]
			text = label_dict[index][2]
			n = 0
			for x in text:
				if x not in Event_List:
 					print "%s is a  new class" % x
					Event_List.append(x)
					f = open(os.path.join(dict_path,x),'w')

				else:
 					f = open(os.path.join(dict_path,x),'a')

				str_convert = wav_dict[index][2][xmin[n]:xmax[n]]
				for i in str_convert:
					print >> f, "%d " % i,
				print >> f
				f.close()
 				n += 1

		wav_list_obj.close()

		events_list_obj = open(os.path.join(dict_path,eventlist),'w')
		Event_List.sort()
		for x in Event_List:
			events_list_obj.write(str(x)+"\n")
		events_list_obj.close()

	print "***************************************"
	print "* Creating List Completed!             *"
	print "***************************************"
if __name__=='__main__':
	if len(sys.argv) != 4:
		exit("Usage: data_list.py <dev_path> <wav_list> <event_list>")
	data_dict(sys.argv[1],sys.argv[2],sys.argv[3])
