'''
This script is intended to calculate the actual Perpendicular Magntetic Field
of the Proton Radiography simulation
'''
import sys
import math
from re import match

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

import numpy as np
