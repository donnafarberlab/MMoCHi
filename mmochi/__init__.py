'''MMoCHi: Multimodal Classifier Hierarchy'''
__version__ = '0.1.2'

from .classifier import terminal_names, predict_proba_cutoff, ground_truth, classifier_setup, classify
from .hierarchy import Hierarchy, gt_defs
from .logger import logg, initiate_log, log_to_file
from .landmarking import landmark_register_adts, density_plot, stacked_density_plots
from .plotting import *

initiate_log()
