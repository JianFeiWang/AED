#!/usr/bin/python

import sys
import os
import re
import numpy as np
import multiprocessing as mp
from sklearn import mixture
#import matplotlib.pyplot as plt
from package.python_speech_features import fbank
#from package.gmm import *

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
	    data = data.strip()
	    data_sum.append(data)
	    nevent += 1
        event_dict[event] = data_sum
        wavdata_obj.close()
        mevents += nevent
    eventlist_obj.close()
    print "The number of events %d" % mevents
    return mevents,event_dict

def medfilt(data):
    data_new = []
    for subdata in data:
        subdata2 = []
        for x in xrange(len(subdata)-2):
     #       print x
    #        raw_input()
            subdata2.append(np.dot(subdata[x:x+3],np.ones(3).T)/3)
        data_new.append(subdata2)
    return data_new

def vad_feature(event_dict,win_ham):
    Vad_vector={}
    for key in event_dict:
	Vad_vector[key] = []

    print 'loading...'

    for key in event_dict:
	ldata = event_dict[key]
	for data in ldata:
   	    adata = []
   	    for sdata in data.split():
                adata.append(int(sdata))
            npoint =  len(adata)
            nframe = (npoint - (32 * 8)) / (10 * 8) + 1
   	    if nframe <= 0:
   	        continue
   	    data_frame = np.zeros((nframe,256))
   	    vector = np.zeros((nframe-1,1))

   	    for n in xrange(nframe):
    	        data_frame[n] = adata[(n * 80) : (256 + n * 80) ]

   	    for n in xrange(nframe-1):
    	        wdata_frame = data_frame[n] * win_ham
   	        E = sum(wdata_frame*wdata_frame)
                vector[n] = [E]

   	    sys.stdout.write('.')
   	    sys.stdout.flush()

            Vad_vector[key].append(vector)

    print "completed!"
    return Vad_vector

def Post_handling(vector):
   #normalization
    means = np.mean(vector,axis = 0)
    stds = np.std(vector,axis = 0)
    vector = (vector - means)/(stds+10E-5)

    return means,stds

def multi_task_list(event_dict,win_ham,ntask):
    Vad_vectors = {}
    ntask = int(ntask)
    key_list = event_dict.keys()
    dict_len = len(key_list)
    names = locals()
    dicts = locals()
    Vdicts = locals()
    len_unit = dict_len / ntask
    for i in xrange(ntask):
	names['key_list%d' % i] = key_list[i*len_unit:(i+1)*len_unit]
	dicts['event_dict%d' % i] = {}
    for key in key_list:
        for i in xrange(ntask):
            if key in names['key_list%d' % i]:
                dicts['event_dict%d' % i][key] = event_dict[key]

    pool = mp.Pool(processes = ntask)
    for i in xrange(ntask):
        Vdicts['Vad_vector%d' % i] = pool.apply_async(vad_feature,(dicts['event_dict%d' % i],win_ham,))
    pool.close()
    pool.join()
    for i in xrange(ntask):
        Vad_vectors.update(Vdicts['Vad_vector%d' % i].get())
    return Vad_vectors

def show_data(arr):
    plt.plot(arr)
    plt.ylabel('some numbers')
    plt.show()

def vadtrain(dict_path,feat_path,tmp_path,eventlist,vadfile,vadfile_norm,vad_statistic,ntask):

    if not os.path.exists(os.path.join(tmp_path,vad_statistic)):
        print "No STATISTIC found!"
        if not os.path.exists(os.path.join(tmp_path,vadfile)):
            print "No VADVECTOR found!"
            Vad_vector = {}
            win_ham = np.hamming(256)
        #   Count = 0
            mevents,event_dict = read_data(dict_path,eventlist)

        ## 32ms hamming 10ms skip
            Vad_vector = multi_task_list(event_dict,win_ham,ntask)
            Vad_vector2 = {}
            for i in event_dict:
                Vad_vector2[i] = []
            vector_v = Vad_vector.values()
            n = 1
            for x in vector_v:
                for d in x:
                    if n == 1:
                        vector = d
                        n = 2
                    else:
                       vector = np.vstack((vector,d))
            means,stds = Post_handling(vector)

            for key,value in Vad_vector.items():
                values = []
                for v in value:
                    va = (v - means)/(stds+10E-05)
                    p = np.ones(1)
                    val = np.dot(va,p)
                    values.append(val)
                    #print values
                    #raw_input()
                Vad_vector2[key].append(medfilt(values))
            f = open(os.path.join(tmp_path,vadfile),"w")
            for key,value in Vad_vector2.items():
                for i in value:
                    for j in i:
                        v = ""
                        for z in j:
                            v = v + " " + str(z)
        	        print >> f, "[%s] %s" % (key,v)
            f.close()
            print "VADVECTOR OK!"
    #   collect the framecount-feature value
 #   normalize
            f = open(os.path.join(tmp_path,vadfile),"r")
            g = open(os.path.join(feat_path,vadfile_norm),"w")
            v = []
            for value in f.readlines():
                value = value.split()
                for value_v in value[1:]:
                    v.append(float(value_v))
            maxvalue = max(v)
            minvalue = min(v)
            n = 3/(maxvalue-minvalue)

            f.seek(0)
            for value in f.readlines():
                value = value.split()
                print >> g, "%s" % value[0],
                for x in value[1:]:
                    k = (float(x) - minvalue) * n
                    print >> g," %.1f" % k,
                print >> g, "\n",

            f.close()
            g.close()

        g = open(os.path.join(feat_path,vadfile_norm),"r")
        k = open(os.path.join(tmp_path,vad_statistic),"w")
        q = np.zeros(31)
        for value in g.readlines():
            value.replace("\n"," ")
            features = re.findall(r' (-?\d+\.?\d*)',value)
            if features:
                for x in features:
                    t = int(float(x) * 10)
                    q[t] += 1
        sumv = sum(q)
        for x in q:
            t = float(x)/sumv *100
            print >> k, "%.2f" % t,
        print >>k, "\n",
        g.close()
        k.close()
        print "STATISTIC OK!"

        print "******************************************"
        print "* hard-decision prepared!                *"
        print "******************************************"

#hard-decision:
#   n-class GMM
    k = open(os.path.join(feat_path,vadfile_norm),"r")
    data = []
    data_label = []
    for value in k.readlines():
        value.replace("\n"," ")
        key = re.findall(r'\[(.+)\]',value)
        features = re.findall(r' (-?\d+\.?\d*)',value)
        if features:
            for x in features:
               data.append([float(x)])
               data_label.append(key[0])
    k.close()

    gmm = mixture.GMM(n_components = 2, covariance_type = 'diag', tol = 0.0001,n_iter = 100, params = "wmc", verbose = 2)
    gmm.fit(data)

    print "******************************************"
    print "* GMM completed!                         *"
    print "******************************************"

    X = data
    RE = gmm.predict(X)
    Score = gmm.score(X)

    f =open("output/Result","w")
    CA1 = 0
    CA2 = 0
    FA  = 0
    MA  = 0
    Z1 = 0
    Z2 = 0
    Z = 0
    for i in xrange(len(RE)):
        print >> f, "%s : %d : %.2f" % (data_label[i],RE[i],Score[i])
        if data_label[i] in ['k1','k2','k3','k4']:
            if RE[i] == 0:
                CA1 += 1
            else:
                MA += 1
            Z1 += 1
        else:
            if RE[i] == 1:
                CA2 += 1
            else:
                FA += 1
            Z2 += 1
        Z += 1
    CA1 = CA1 / float(Z1) *100
    MA = MA / float(Z1) *100
    FA = FA / float(Z) *100
    print >> f, "FA = %.2f\nMA = %.2f\nCA = %.2f" % (FA,MA,CA1)
    f.close()


    return gmm

def vad_frame(gmm,data):


#    win_ham = np.hamming(256)
#    adata = []
#    for sdata in testdata.split():
#        adata.append(int(sdata))
#    npoint =  len(adata)
#    nframe = (npoint - (32 * 8)) / (10 * 8) + 1
#    data_frame = np.zeros((nframe,256))
#    vector = np.zeros((nframe-1,5))
#
#    for n in xrange(nframe):
#        data_frame[n] = adata[(n * 80) : (256 + n * 80) ]
#
#    for n in xrange(nframe-1):
#        wdata_frame = data_frame[n] * win_ham
#        wdata_frame2 = data_frame[n+1] * win_ham
#        h,c,Gp = time_feature_extract(wdata_frame,win_ham)
#        Phps,SFp = frequency_feature_extract(wdata_frame,wdata_frame2)
#        vector[n] = [h,c,Gp,Phps,SFp]
#
#    values = []
#    for v in vector:
#        va = (v - means)/(stds+10E-05)
#        val = np.dot(va,M.T)
#        values.append(val)
#
#    k =(value - minvalue)*n
#    q = np.zeros(31)
#    for w in k:
#        t = int(float(x)*10)
#        q[t]+=1
#    sumv = sum(q)
#    data = q/sumv*100
#    Re = gmm.predit(data)
#    Score = gmm.score(data)

    return Re,Score

if __name__=="__main__":
    '''
    vad test part need do more things
    '''


    if len(sys.argv) != 9:
        exit("Usage: vad.py <dict_path> <feature_path> <tmp_path> <eventlist> <vadvector> <vadvector_norm> <vad_feature> <ntask>")
    gmm = vadtrain(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])
