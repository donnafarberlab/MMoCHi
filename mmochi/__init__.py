'''MMoCHi: Multimodal Classifier Hierarchy'''
__version__ = '0.1.4'

from .classifier import terminal_names, predict_proba_cutoff, ground_truth, classifier_setup, classify, idenitfy_group_markers
from .hierarchy import Hierarchy, gt_defs
from .logger import logg, initiate_log, log_to_file
from .landmarking import landmark_register_adts, density_plot, stacked_density_plots, density_plot_total, \
                         save_peak_overrides, load_peak_overrides, update_peak_overrides, update_landmark_register
from .plotting import plot_tree, feature_importances, plot_confusion, plot_confidence, plot_important_features

from . import classifier
from . import hierarchy
from . import landmarking
from . import logger
from . import plotting
from . import thresholding
from . import utils

initiate_log()