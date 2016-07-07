# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:35:58 2016

@author: WolfeTM
"""
import numpy as np

b1 = np.zeros((4,5))
b2 = np.ones((4,3))

pos_v, pos_h = 2, 3  # offset

v_range1 = slice(max(0, pos_v), min(pos_v + b2.shape[0], b1.shape[0]))
h_range1 = slice(max(0, pos_h), min(pos_h + b2.shape[1], b1.shape[1]))

v_range2 = slice(max(0, -pos_v), min(-pos_v + b1.shape[0], b2.shape[0]))
h_range2 = slice(max(0, -pos_h), min(-pos_h + b1.shape[1], b2.shape[1]))

b1[v_range1, h_range1] += b2[v_range2, h_range2]