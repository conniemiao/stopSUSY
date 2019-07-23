#!/usr/bin/env python

# Train MLP Classifier on some signal and bkgd data, fit it, test it on some
# other signal and bkgd data, and plot the ROC curve.

print "Importing modules."
# from sklearn.neural_network import MLPClassifier
# from sklearn import tree
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import normalize
from sklearn import datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print "Importing training data."
iris = datasets.load_iris()
XTrain = np.concatenate((iris.data[:25,:2], iris.data[50:75,:2], iris.data[100:125,:2]))
yTrain = np.concatenate((iris.target[:25], iris.target[50:75], iris.target[100:125]))
# XTrain = np.concatenate((iris.data[:1,:2], iris.data[55:56,:2]))
# yTrain = np.concatenate((iris.target[:1], iris.target[55:56]))
print XTrain
print yTrain
for i, val in enumerate(yTrain):
    if val > 0: yTrain[i] = 1

print "Training."
# clf = MLPClassifier(hidden_layer_sizes=(2,2))
# clf = tree.DecisionTreeClassifier()
clf = GradientBoostingClassifier()
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
print pd.DataFrame([clf.feature_importances_], columns=iris.feature_names[:2])
plt.figure()
plt.plot(fpr, tpr, label = 'auc=%0.2f' % roc_auc)
plt.title("ROC curve GBC iris")
plt.legend(loc = 'lower right')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
axes = plt.gca()
axes.set_xlim([0,1])
axes.set_ylim([0,1])
plt.show()
