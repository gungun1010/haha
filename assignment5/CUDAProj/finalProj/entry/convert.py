#!/usr/bin/python
import commands
import time 
import struct
import os
import array
import binascii

#number of meaningful bytes for GPU
PATTERN_BYTES = 10

def parse(ndb, gpusig, gpuvirus):
    #sigs = []
    #virus = []
    f = open(ndb)
    fs = open(gpusig,'wb')
    fv = open(gpuvirus,'w')
    lines = f.readlines()
    f.close()
    
    offset = PATTERN_BYTES*2

    for idx in range(len(lines)):
        #print idx
        info = lines[idx].split(":")
        #info[0] is virus's name
        #info[3] is virus signature
    
        sigs = (info[3][:offset])
        virus = (info[0])
        #print sigs
        try:
            sigsBytes = binascii.a2b_hex(sigs)
            fs.write (sigsBytes)
            fv.write (virus+os.linesep)
        except TypeError:
           next 
        #print sigs[idx]
        #print virus[idx]
        #time.sleep(1)
    fs.close()
    fv.close()


def main():
    print "start converting for main"
    ndb = "/home/leon/clamav/share/clamav/mainPack/main.ndb"
    gpusig = "/home/leon/clamav/share/clamav/mainPack/mainGPUsig.bin"
    gpuvirus = "/home/leon/clamav/share/clamav/mainPack/mainGPUvirus.ndb"

    parse(ndb, gpusig, gpuvirus)
    
    print "start converting for daily"
    ndb = "/home/leon/clamav/share/clamav/dailyPack/daily.ndb"
    gpusig = "/home/leon/clamav/share/clamav/dailyPack/dailyGPUsig.bin"
    gpuvirus = "/home/leon/clamav/share/clamav/dailyPack/dailyGPUvirus.ndb"
    parse(ndb, gpusig, gpuvirus)
    
main()

