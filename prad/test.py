# -*- coding: utf-8 -*-⏊

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

x = np.array(range(1,11))
y = np.array(range(20,26))
X, Y = np.meshgrid(x,y)
print X, X.shape
print Y, Y.shape
