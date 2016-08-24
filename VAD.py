#!/usr/bin/python

##############################################
#  LIBRARY AND FUNCTIONS                     #
##############################################
import os
import sys
import re
import datetime
import numpy as np
import wave as wv
from data_vad import multi_task_list as mtl
from data_vad import medfilt as mf
from sklearn import mixture as mi
from logger import *

#############################################
#  DEFINED VARIABLES                        #
#############################################
CV_ADDRESS = 'data/cv/'
PARAMS_ADDRESS = 'params/'
TMP_ADDRESS = 'tmp/'
FEAT_ADDRESS = 'feature/'
OUTPUT_ADDRESS = 'output/'
LOG_ADDRESS =  'log/'
VADWAV_ADDRESS = 'data/vw'
VAD_FEATURE = 'vadvector_norm'
VAD_RESULT = 'vad_result'
VAD_SCORE_THRED = 0.3
VAD_R_TRUE = 0
VAD_R_FALSE = 1
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


#  OPEN THE LOGFILE
#now = datetime.datetime.now()
#DATE = now.strftime('%Y-%m-%d')
#TIME = now.strftime('%H:%M:%S')
#flog = open(os.path.join(LOG_ADDRESS,DATE),'w+')
#print >> flog, "[%s]" % TIME

#  READ THE VAD MODEL AND PARAMS
#print >> flog, "READ VAD MODEL PARAMS.............",
log = LOGGER()
log.UPLoad("READ VAD MODEL PARAMS.............",)
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

log.UPLoad("SUCCESS")
log.PUSH()
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
log = LOGGER()
log.UPLoad("TRAIN THE GMM MODEL............",)
#print >> flog, "TRAIN THE GMM MODEL............",
fvadp = open(os.path.join(PARAMS_ADDRESS,'par_vad_gmm'),'r')
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

gmm = mi.GMM(n_components = 2, covariance_type = 'diag', tol = 0.0001,   n_iter = 100, params = "wmc", verbose = 0)
gmm.fit(data)
log.UPLoad("SUCCESS")
log.PUSH()


#  READ THE CVDATA AND EXTRACT THE VADFEATURE

fcv_list = open(os.path.join(CV_ADDRESS,'cv.list'),"r")
for cv in fcv_list.readlines():
#    print >> flog, "%s DETECTIMG..........." % cv

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


    log = LOGGER()
    log.UPLoad("VADFEATURE EXACTING...........",)
    VAD_vector = mtl(DATADICT,HAMMING,nkey)   # THE LAST FRAME IS LOST HERE, SO WHEN DO THIS PART, WE OLEY GET THE FRAME -1 BACK
    plog.UPLoad("SUCCESS")
    log.PUSH()
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
    log = LOGGER()
    log.UPLoad("PCA...........",)
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

    log.UPLoad("SUCCESS")
    log.PUSH()
#3        print vb
#        raw_input()
#        VAD_vector_norm = VAD_vector
#    fcvv = open(os.path.join(TMP_ADDRE))



#   VAD TASK
    log = LOGGER()
    log.UPLoad("VAD...........",)
    Label = {}
    VADResult ={}
    for key in VAD_vector3.keys():
        Label[key] = [[],[]]
        VADR[key] = []
      #  print VAD_vector3[key]
       # raw_input()
    for key,value in VAD_vector3.items():
  #      log.UPLoad('value0',type(value),len(value))
        files = []
        for v in value:
#            log.UPLoad('value1',type(v),len(value))
            for vx in v:
 #               log.UPLoad('value2',type(vx),len(value))
                for vxx in vx:
                    X = vxx
                    RE = gmm.predict(X)
                    Score = gmm.score(X)
                    if (Score < VAD_SCORE_THRED) and (RE == VAD_R_FALSE):
                       REs = VAD_R_TRUE
                    else:
                       REs = RE
                    Label[key][0].append(RE)
                    Label[key][1].append(RE)
                    files.append(RE)
        VADR[key].append(files)

    fvadr = open(os.path.join(OUTPUT_ADDRESS,VAD_RESULT+cv),"w")

    for key,value in Label.items():
        Len = len(value[0])
        for l in xrange(Len):
            print >> fvadr, "[%s] %d  %f" % (key,value[0][l],value[1][l])

    log.UPLoad("SUCCESS")
    log.PUSH()

    fvadr.close()

    for key,files in DATADICT.items():
        for fs in files:
            nframe = len(fs) - 3
            wav = []
            for n in xrange(nframe):
                wav.append(double(VADR[key][n]) * double(fs[n]))

            nZeros = 0
            BWAVOPEN = 0
            BWAVW = 0
            fwav = []
            ff = []
            for w in wav:
                if(w != 0):
                    if(BWAVOPEN == 0):
                        BWAVOPEN = 1
                    BWAVW = 1
                    nZeros = 0
                if(w == 0):
                    nZeros += 1
                    if(nZeros == 10):
                        BWAVW = 0;
                if(BWAVW):
                    if(BWAVOPEN == 0):
                        fwav.append(ff)
                        ff = []
                    ff.append(w)
                else:
                    if(BWAVOPEN):
                        BWAVOPEN = 0;

            num = 0;
            for fw in fwav:
                L = len(fw)
                f= wv.open(os.path.join(VADWAV_ADDRESS, key + str(num)),'w')
                num += 1;
                f.setparams((1,2,8000,L,NONE,key))
                f.writeframesraw(fw)
                f.close()

    fcv.close()
    print "OK"

fcv_list.close()

