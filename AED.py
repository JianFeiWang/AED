#!/usr/bin/python

##############################################
#  LIBRARY AND FUNCTIONS                     #
##############################################
import os
import sys
import re
import numpy as np
from data_vad import multi_task_list as mtl
from data_vad import medfilt as mf
from sklearn import mixture as mi

#############################################
#  DEFINED VARIABLES                        #
#############################################
CV_ADDRESS = 'data/cv/'
PARAMS_ADDRESS = 'params/'
TMP_ADDRESS = 'tmp/'
FEAT_ADDRESS = 'feature/'
OUTPUT_ADDRESS = 'output/'
VAD_FEATURE = 'vadvector_norm'
VAD_RESULT = 'vad_result'
HAMMING = np.hamming(256)


#############################################
#  IMPORTANT COMMAN VARIABLES               #
#############################################
KEY = ''
DATADICT = {}
VADFEATUR = {}
VAD_PARAMS = []
NFRAME = 0
DATA_FRAME_MAP = {}

############################################
#  MAIN FUNCTION                           #
############################################

#  READ THE VAD MODEL AND PARAMS
fvadp = open(os.path.join(PARAMS_ADDRESS,'par_vad_gmm'),'r')
for f in fvadp.readlines():
    f = f.strip()
    f = re.findall(r'(\-?\d+\.?\d*e?\+?\-?\d*)',f)
    flen = len(f)
    if flen > 1:
        n = np.zeros(flen)
        for x in xrange(flen):
            n[x] = float(f[x])
    else:
        n = float(f[0])
    VAD_PARAMS.append(n)

vad_params_mul = VAD_PARAMS[0]
vad_params_bias = VAD_PARAMS[1]
vad_params_M = VAD_PARAMS[2]
vad_params_means = VAD_PARAMS[3]
vad_params_stds = VAD_PARAMS[4]

print "VAD PARAMATERS:"
print "MEANS:",
print vad_params_means
print "STDS:",
print vad_params_stds
print "MODMATRIX:",
print vad_params_M
print "MULTI:",
print vad_params_mul
print "BIAS:",
print vad_params_bias
print "\n",

#   TRAIN THE GMM MODEL
k = open(os.path.join(FEAT_ADDRESS,VAD_FEATURE),"r")
data = []
data_label = []
for value in k.readlines():
    value.replace("\n"," ")
    key = re.findall(r'\[(.+)\]',value)
    features = re.findall(r' (-?\d+\.?\d*)',value)
    if features:
        for x in features:
           data.append([float(x)])
     #      print data
     #      raw_input()
           data_label.append(key[0])
k.close()

gmm = mi.GMM(n_components = 2, covariance_type = 'diag', thresh = 1e-3, min_covar = 1e-7, n_iter = 100, params = "wmc")
gmm.fit(data)

#  READ THE CVDATA AND EXTRACT THE VADFEATURE
fcv_list = open(os.path.join(CV_ADDRESS,'cv.list'),"r")

for cv in fcv_list.readlines():

    cv = cv.strip()
    print cv
    print "***********"
    fcv = open(os.path.join(CV_ADDRESS,cv),'r')
    nkey = 0

    n = 0
    DATADICT = {}
    for data in fcv.readlines():

        data = data.strip()
        key = re.findall(r'\[\'([a-z]+\d?)\'\]',data)
        if key:
            KEY = key[0]
            if not KEY in DATADICT.keys():
                nkey += 1
                print "%s\n" % KEY,
                DATADICT[KEY] = []
            continue
#            print "..."
        DATADICT[KEY].append(data[2:-2])
        n += 1
    print "The number of events : %d" % n

    VAD_vector = mtl(DATADICT,HAMMING,nkey)   # THE LAST FRAME IS LOST HERE, SO WHEN DO THIS PART, WE OLEY GET THE FRAME -1 BACK

#1    print VAD_vector[sli]
#    raw_input()
 #   fvv = open(os.path.join(TMP_ADDRESS,'vadvectorT'),'w')
 #   for key,value in VAD_vector.items():
#        print key
#        raw_input()
#
#        print value
#        raw_input()
#
 #       for v in value:
 #           print >> fvv, "[%s]" % key,
 #           for x in v:
 #               print x
 #               print >> fvv, " %f" % float(x),
 #           print >> fvv, "\n",
 #   fvv.close()

    VAD_vector2 = {}
    VAD_vector3 = {}

    for key in VAD_vector.keys():
        VAD_vector2[key] = []
        VAD_vector3[key] = []

#    NFRAME = len(VAD_vector)
#    for n in xrange(NFRAME):
#        DATA_FRAME_MAP[n] = DATA

    for key,value in VAD_vector.items():
        vc = []
        for v in value:
    #        print v
     #       raw_input()
            va = (v - vad_params_means)/(vad_params_stds+10E-05)
            vb = np.dot(va,vad_params_M.T)
            vc.append(vb)
        vd = mf(vc)  #THE LAST TWO FRAME ARE LOST HERE, SO WHEN DO THIS PART, WE ONLY GET THE FRAME -2 BACK
        VAD_vector2[key].append(vd)
#2        print vd
#        raw_input()


    for key,value in VAD_vector2.items():
        vb = []
        for v in value:
           # print len(value)
#            raw_input()
            for vs in v:
                vx = np.zeros(len(vs))
                for i in xrange(len(vs)):
                    vx[i] = float(vs[i])
              #  print vx
              #  print type(vx)
              #  raw_input()
                va = vx * vad_params_mul - vad_params_bias
                vb.append(va)
               # print vb
               # print len(vb)
               # raw_input()
        VAD_vector3[key].append(vb)
#3        print vb
#        raw_input()
#        VAD_vector_norm = VAD_vector
#    fcvv = open(os.path.join(TMP_ADDRE))



#   VAD TASK
    Label = {}
    for key in VAD_vector3.keys():
        Label[key] = []
      #  print VAD_vector3[key]
       # raw_input()
    for key,value in VAD_vector3.items():
        kv = np.zeros(len(value))
        for v in value:
            for vx in v:
                for vxx in vx:
                    X = vxx
                #    print X
                #    raw_input()
                    RE = gmm.predict(X)
                    Score = gmm.score(X)
                    Label[key].append(RE)

    fvadr = open(os.path.join(OUTPUT_ADDRESS,VAD_RESULT+cv),"w")

    for key,value in Label.items():
        for v in value:
            print >> fvadr, "[%s]" % key,
            for vx in v:
                print >>fvadr, " %d" % vx,
            print >>fvadr, "\n",

    fvadr.close()

    fcv.close()
    print "OK"


fcv_list.close()


