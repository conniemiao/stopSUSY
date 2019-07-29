#!/usr/bin/env python

# Train a machine learning algorithm on some signal and bkgd data, fit it, 
# test it on some other signal and bkgd data, and plot the ROC curve.

print "Importing modules."
# from sklearn.neural_network import MLPClassifier
# from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import normalize
from sklearn.utils import shuffle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

branches = ['met_pt','met_phi', 'njets','nbtag','jet_ht', 'lep1_pt', 'lep2_pt', \
        'lep1_mt', 'lep2_mt', 'mt_tot', 'mt_sum', 'm_eff']

print "Importing training data."
sigTrainAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/sig0_elel.hdf"
bkgdTrainAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/stopCut_all_Bkgd_TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf"
sigTrainDF = pd.read_hdf(sigTrainAdr, "tSig0")
bkgdTrainDF = pd.read_hdf(bkgdTrainAdr, "tBkgd")
print "num sig evts:", sigTrainDF.shape[0], "num bkgd evts tot:", bkgdTrainDF.shape[0]
bkgdTrainDF = bkgdTrainDF.iloc[:100*sigTrainDF.shape[0],:]
# bkgdTrainDF = pd.read_hdf(bkgdTrainAdr, "tBkgd")
XTrain = pd.concat([bkgdTrainDF, sigTrainDF])
# XTrain = normalize(XTrain, norm='max', axis=0)
yTrain = [1]*sigTrainDF.shape[0]+[0]*bkgdTrainDF.shape[0]
# print XTrain.head(15)
# print yTrain[:15]
# print "shuffle"
XTrain, yTrain = shuffle(XTrain, yTrain)
# print XTrain.head(15)
# print yTrain[:15]

print "Training."
# clf = MLPClassifier(hidden_layer_sizes=(12,2))
# clf = DecisionTreeClassifier()
clf = GradientBoostingClassifier()
clf.fit(XTrain, yTrain)


importancesDF = pd.DataFrame([clf.feature_importances_], columns=branches)
print "Feature importances:"
print importancesDF

print "Importing testing data."
# sigTestAdr = sigTrainAdr
# bkgdTestAdr = bkgdTrainAdr
# sigTestDF = pd.read_hdf(sigTestAdr, "tSig0")
# bkgdTestDF = pd.read_hdf(bkgdTestAdr, "tBkgd")

sigTestAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/sig1_elel.hdf"
bkgdTestAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/HLLHC/pandas_files/stopCut_all_Bkgd_TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf"
sigTestDF = pd.read_hdf(sigTestAdr, "tSig1")
bkgdTestDF = pd.read_hdf(bkgdTestAdr, "tBkgd").iloc[101*sigTrainDF.shape[0]:102*sigTrainDF.shape[0],:]
# bkgdTestDF = pd.read_hdf(bkgdTestAdr, "tBkgd")

XTest = pd.concat([bkgdTestDF, sigTestDF])
# XTest = normalize(XTest, norm='max', axis=0)
yTest = [1]*sigTestDF.shape[0]+[0]*bkgdTestDF.shape[0]

print "Testing."
print "Got", clf.predict_proba(XTest)[:300,1]
print "Should be", yTest[:300]
fpr, tpr, thresholds = roc_curve(yTest, clf.predict_proba(XTest)[:,1])
# area under curve
roc_auc = auc(fpr, tpr)

print "Plotting."
importancesDF = pd.DataFrame([clf.feature_importances_], columns=branches)
print "Feature importances:"
print importancesDF
plt.figure()
plt.plot(fpr, tpr, label = 'auc=%0.2f' % roc_auc)
plt.title("ROC curve GBC ttbar bkgd vs. sig")
plt.legend(loc = 'lower right')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
axes = plt.gca()
axes.set_xlim([0,1])
axes.set_ylim([0,1])
plt.show()
