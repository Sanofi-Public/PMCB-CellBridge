import matplotlib
matplotlib.use('agg')
import scrublet as scr
import scipy.io
import numpy as np
import os
from scipy.sparse import csc_matrix

def scrublet_py(i, j, val, dim):
	data = csc_matrix((val, (i, j)), shape = dim)
	scrub = scr.Scrublet(data, random_state=123)
	doublet_scores, predicted_doublets = scrub.scrub_doublets()
	return(doublet_scores, scrub.doublet_scores_sim_) # predicted_doublets, 

