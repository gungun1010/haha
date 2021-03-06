#!/usr/bin/python

import os
import h5py
import sys
import shutil
import sklearn
import tempfile
import numpy as np
import pandas as pd
import sklearn.datasets
import sklearn.linear_model
import matplotlib.pyplot as plt
import subprocess
from subprocess import call
import commands
import time 
import xxhash
import struct

NUM_COLS = 60
NUM_ITR = 100

class NestedDict(dict):
    def __getitem__(self, key):
        if key in self: return self.get(key)
        return self.setdefault(key, NestedDict())

def split_hex(text):
    text = text.zfill(len(text) + len(text) % 16)
    return ' '.join(text[i: i+16] for i in range(0, len(text), 16))

def sklearnSGD():
    global trainData
    global trainLabel 
    global testData
    global testLabel
    
    trainData, trainLabel = sklearn.datasets.make_classification(
        n_samples=14300, n_features=NUM_COLS, n_redundant=0, n_informative=2, 
        n_clusters_per_class=2, hypercube=False, random_state=0
    )

    # Split into train and test
    trainData, testData, trainLabel, testLabel = sklearn.cross_validation.train_test_split(trainData, trainLabel)

    # Visualize sample of the data
    #ind = np.random.permutation(trainData.shape[0])[:1000]
    #df = pd.DataFrame(trainData[ind])
    #_=pd.scatter_matrix(df, figsize=(9, 9), diagonal='kde', marker='o', s=40, alpha=.4, c=y[ind])

    # Train and test the scikit-learn SGD logistic regression.
    #clf = sklearn.linear_model.SGDClassifier(
     #   loss='log', n_iter=1000, penalty='l2', alpha=1e-3, class_weight='auto')

    #clf.fit(trainData, trainLabel)
    #testLabel_pred = clf.predict(testData)
    #print('Accuracy: {:.3f}'.format(sklearn.metrics.accuracy_score(testLabel, testLabel_pred)))

def parseProcess():
    global trainData, trainLabel, testData, testLabel
    
    trainData = np.array([],dtype=float)
    procInfoHex = []
    #trainLabel = np.array([],dtype=int)
    testData = np.array([], dtype=float)
    #testLabel = np.array([], dtype=int)

    #entry[0] = ['USER', 'PID', '%CPU', '%MEM', 'VSZ', 'RSS', 'TTY', 'STAT', 'START', 'TIME', 'COMMAND']
    #USER(0), %CPU(2), %MEM(3), VSZ(4), RSS(5), TIME(9), CMD(10) are important attribute
    init = 1
    iteration = 0

    while iteration < NUM_ITR: #100 iteration, about 2000 samples
        shDump = commands.getoutput("ps aux --sort=-pcpu") # list of all processes in stream
        entry = shDump.split('\n') #split the stream into arraies
        #print len(entry)

        for indx in range(len(entry)):
            hexList = []
            if indx != 0:
                procInfo = entry[indx].split() #split array into elements
                
                #remove meaningless contents
                procInfo.pop(1);
                procInfo.pop(5);
                procInfo.pop(6);
                #procInfo = ['USER', '%CPU', '%MEM', 'VSZ', 'RSS','STAT','TIME', 'COMMAND','param0', 'param1', 'param2', 'param3', 'param4', .... 'paramN']
                #convert plain text string to 8-byte element a pop
                for s in procInfo:
                    #if an attribute is more than 8 bytes, we split it
                    procInfoHex = split_hex(s.encode('hex'))
                    procInfoHex = procInfoHex.split()
                    hexList= hexList + procInfoHex
                
                
                #the buffer is 60 doubles wide, so padding the remainders with ' '
                while len(hexList) < NUM_COLS:
                    hexList.append(' '.encode('hex'))
                
                #covert Hex to double
                procInfoHash = [ struct.unpack('<d', s.zfill(16).decode('hex'))[0] for s in hexList]

                if indx == 1 and init == 1:
                    init=0 #only init numpy once
                    
                    #procInfo[7] is the command, use it as label
                    trainData=np.hstack((trainData, np.float_(procInfoHash)));
                    #trainLabel=np.append(trainLabel, np.int_(procInfoHash[7]));

                else:
                    trainData=np.vstack((trainData, np.float_(procInfoHash)));
                    #trainLabel=np.append(trainLabel, np.int_(procInfoHash[7]));
        iteration += 1
    
    # Split into train and test
    trainData, testData= sklearn.cross_validation.train_test_split(trainData)

def genH5():
    global dirname, train_filename, test_filename

    # Write out the data to HDF5 files in a temp directory.
    # This file is assumed to be caffe_root/examples/hdf5_classification.ipynb
    dirname = os.path.abspath('/home/leon/enas907/proj/data')
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    train_filename = os.path.join(dirname, 'train.h5')
    test_filename = os.path.join(dirname, 'test.h5')

    # HDF5DataLayer source should be a file containing a list of HDF5 filenames.
    # To show this off, we'll list the same data file twice.
    with h5py.File(train_filename, 'w') as f:
        f['data'] = trainData
        f['label'] = trainLabel.astype(np.float32)
    with open(os.path.join(dirname, 'train.txt'), 'w') as f:
        f.write(train_filename + '\n')
        f.write(train_filename + '\n')

def zipH5():
    # HDF5 is pretty efficient, but can be further compressed.
    global test_filename

    comp_kwargs = {'compression': 'gzip', 'compression_opts': 1}
    with h5py.File(test_filename, 'w') as f:
        f.create_dataset('data', data=testData, **comp_kwargs)
        f.create_dataset('label', data=testLabel.astype(np.float32), **comp_kwargs)
    with open(os.path.join(dirname, 'test.txt'), 'w') as f:
        f.write(test_filename + '\n')

def trainDNN():
    call(["caffe", "train","-solver","/home/leon/enas907/proj/solver.prototxt"])

def main():
    global trainData, trainLabel, testData, testLabel, dirname, train_filename, test_filename 
    
    sklearnSGD()
    print "this is working dataset -----------------------------------------"
    print "trainData: "
    print trainData
    print trainData.shape
    print "trainLabel: "
    print trainLabel
    print trainLabel.shape
    print "testData: "
    print testData
    print testData.shape
    print "testLabel: "
    print testLabel
    print testLabel.shape

    parseProcess()
    print "this is working dataset -----------------------------------------"
    print "trainData: "
    print trainData
    print trainData.shape
    print "trainLabel: "
    print trainLabel
    print trainLabel.shape
    print "testData: "
    print testData
    print testData.shape
    print "testLabel: "
    print testLabel
    print testLabel.shape
    #sys.exit()


    genH5()

    zipH5()
    
    #FIXME try train tomorrow
    trainDNN()

main()
