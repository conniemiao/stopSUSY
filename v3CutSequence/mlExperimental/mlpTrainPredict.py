# Train MLP Classifier on some signal and bkgd data, fit it, test it on some
# other signal and bkgd data, and plot the ROC curve.

print "Importing"
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import normalize
import pandas as pd
import matplotlib.pyplot as plt

print "Training"
sigTrainAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/pandas_files/sig0_elel.hdf"
bkgdTrainAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/pandas_files/all_Bkgd_DY0Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf"
sigTrainDF = pd.read_hdf(sigTrainAdr, "tSig0")
bkgdTrainDF = pd.read_hdf(bkgdTrainAdr, "tBkgd")
XTrain = normalize(pd.concat([bkgdTrainDF, sigTrainDF]), norm='max', axis=0)
yTrain = [1]*sigTrainDF.shape[0]+[0]*bkgdTrainDF.shape[0]

clf = MLPClassifier(hidden_layer_sizes=(4,2))
clf.fit(XTrain, yTrain)

print "Testing"
sigTestAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/pandas_files/sig0_elel.hdf"
bkgdTestAdr = "/afs/cern.ch/work/c/cmiao/private/myDataSusy/pandas_files/all_Bkgd_DY0Jets_MLL-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_elel.hdf"
sigTestDF = pd.read_hdf(sigTestAdr, "tSig0")
bkgdTestDF = pd.read_hdf(bkgdTestAdr, "tBkgd")
XTest = normalize(pd.concat([bkgdTestDF, sigTestDF]), norm='max', axis=0)
# XTest = pd.concat([bkgdTestDF, sigTestDF])
yTest = [1]*sigTestDF.shape[0]+[0]*bkgdTestDF.shape[0]

fpr, tpr, thresholds = roc_curve(yTest, clf.predict_proba(XTest)[:,1])
# area under curve
roc_auc = auc(fpr, tpr)

print "Plotting"
plt.figure()
plt.plot(fpr, tpr, label = 'auc=%0.2f' % roc_auc)
plt.title("ROC curve")
plt.legend(loc = 'lower right')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()
