'''MMoCHi: Multimodal Classifier Hierarchy'''
__version__ = '0.2.1'

from .classifier import terminal_names, hc_threshold, classifier_setup, classify, identify_group_markers, DEBUG_ERRORS
from .hierarchy import Hierarchy, hc_defs, Classification, Subset
from .logger import logg, _initiate_log, log_to_file
from .landmarking import landmark_register_adts, density_plot, stacked_density_plots, density_plot_total, \
                         save_peak_overrides, load_peak_overrides, update_peak_overrides, update_landmark_register
from .plotting import plot_tree, feature_importances, plot_confusion, plot_confidence, plot_important_features

from .utils import DATA_KEY, BATCH_KEY, preprocess_adatas, umap_thresh, umap_interrogate_level

from . import classifier
from . import hierarchy
from . import landmarking
from . import logger
from . import plotting
from . import thresholding
from . import utils

_initiate_log()