#!/usr/bin/env python

# Train MLP Classifier on some signal and bkgd data, fit it, test it on some
# other signal and bkgd data, and plot the ROC curve.

print "Importing modules."
# from sklearn.neural_network import MLPClassifier
from sklearn import tree
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import normalize
from sklearn import datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print "Importing training data."
iris = datasets.load_iris()
print iris.feature_names
XTrain = np.concatenate((iris.data[:25,:2], iris.data[50:75,:2], iris.data[100:125,:2]))
yTrain = np.concatenate((iris.target[:25], iris.target[50:75], iris.target[100:125]))
# XTrain = np.concatenate((iris.data[:1,:], iris.data[55:56,:]))
print XTrain
# yTrain = np.concatenate((iris.target[:1], iris.target[55:56]))
print yTrain
for i, val in enumerate(yTrain):
    if val > 0: yTrain[i] = 1

print "Training."
# clf = MLPClassifier(hidden_layer_sizes=(2,2))
clf = tree.DecisionTreeClassifier()
# clf.fit(XTrain, yTrain)
clf.fit(XTrain, yTrain)

print "Importing testing data."
XTest = np.concatenate((iris.data[25:50,:2], iris.data[75:100,:2], iris.data[125:150,:2]))
yTest = np.concatenate((iris.target[25:50], iris.target[75:100], iris.target[125:150]))
print yTest
for i, val in enumerate(yTest):
    if val > 0: yTest[i] = 1

print "Testing."
print clf.predict_proba(XTest)[:,1]
fpr, tpr, thresholds = roc_curve(yTest, clf.predict_proba(XTest)[:,1])
# area under curve
roc_auc = auc(fpr, tpr)

print "Plotting."
print clf.feature_importances_
plt.figure()
plt.plot(fpr, tpr, label = 'auc=%0.2f' % roc_auc)
plt.title("ROC curve")
plt.legend(loc = 'lower right')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()
