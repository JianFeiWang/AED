#!usr\bin\python
import os
import datetime

LOGADDRESS = 'log/'

class LOGGER:
    def GetDT(self):
        now = datetime.datetime.now()
        self.DATE = now.strftime('%Y-%m-%d')
        self.TIME = now.strftime('%H:%M:%S')
    def __init__(self,addr = 'log/'):
        self.GetDT()
        self.ADDRESS = addr
        self.flog = open(os.path.join(self.ADDRESS,self.DATE),'a')

    def UPLoad(self,*Info):
        self.GetDT()
        print >> self.flog, "[%s]" % self.TIME,
        for inf in Info:
            print >> self.flog, "%r " % inf,
        print >> self.flog, '\n'

    def PUSH(self):
        self.flog.close()
