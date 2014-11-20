#!/usr/bin/python

#%matplotlib inline
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

class NestedDict(dict):
    def __getitem__(self, key):
        if key in self: return self.get(key)
        return self.setdefault(key, NestedDict())

def sklearnSGD():
    global trainData
    global trainLabel 
    global testData
    global testLabel
    
    trainData, trainLabel = sklearn.datasets.make_classification(
        n_samples=10000, n_features=4, n_redundant=0, n_informative=2, 
        n_clusters_per_class=2, hypercube=False, random_state=0
    )

    # Split into train and test
    trainData, testData, trainLabel, testLabel = sklearn.cross_validation.train_test_split(trainData, trainLabel)

    print trainData.shape
    print testData.shape
    print trainLabel.shape
    print testLabel.shape
    sys.exit()
    # Visualize sample of the data
    #ind = np.random.permutation(trainData.shape[0])[:1000]
    #df = pd.DataFrame(trainData[ind])
    #_=pd.scatter_matrix(df, figsize=(9, 9), diagonal='kde', marker='o', s=40, alpha=.4, c=y[ind])

    # Train and test the scikit-learn SGD logistic regression.
    clf = sklearn.linear_model.SGDClassifier(
        loss='log', n_iter=1000, penalty='l2', alpha=1e-3, class_weight='auto')

    clf.fit(trainData, trainLabel)
    testLabel_pred = clf.predict(testData)
    print('Accuracy: {:.3f}'.format(sklearn.metrics.accuracy_score(testLabel, testLabel_pred)))

def parseProcess():
    global trainData, trainLabel, testData, testLabel
    
    #entry[0] = ['USER', 'PID', '%CPU', '%MEM', 'VSZ', 'RSS', 'TTY', 'STAT', 'START', 'TIME', 'COMMAND']
    #USER(0), %CPU(2), %MEM(3), VSZ(4), RSS(5), TIME(9), CMD(10) are important attribute
    shDump = commands.getoutput("ps aux --sort=-pcpu") # list of all processes in stream
    #print shDump 
    entry = shDump.split('\n') #split the stream into arraies
    print len(entry)
    vals = entry[1].split() #split array into elements

    #if there is total number of parameters is less than 5, fill in with blanks
    while len(vals) <= 15:
        vals.append(' ')

    vals.pop(1);
    vals.pop(5);
    vals.pop(6);
    #vals = ['USER', '%CPU', '%MEM', 'VSZ', 'RSS','STAT','TIME', 'COMMAND','param0', 'param1', 'param2', 'param3', 'param4']

    trainData=np.array(vals);
    print trainData

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
    
    parseProcess()
    sys.exit()
    sklearnSGD()

    genH5()

    zipH5()

    trainDNN()

    print "end of program\n"

main()
