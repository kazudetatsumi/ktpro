#!/usr/bin/env python
import numpy as np
from sklearn.model_selection import KFold 

X = np.array([["DIAMOND"], ["SNEAK"], ["A"], ["BOY"], ["MEETS"], ["GIRL"], ["IN"], ["GARDEN"], ["TREE"]])
#y = np.array([1, 2, 3, 4])
kf = KFold(n_splits=3, shuffle=True)
#for train_index, test_index in kf.split(X):
#       print("TRAIN:", X[train_index], "TEST:", X[test_index])
print(kf.split(X))

