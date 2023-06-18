import numpy as np
import pandas as pd
import itertools
import treelib
import pickle
import gc
import re
import sklearn.calibration
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable, List, Set, Dict
import anndata

from . import utils
from .logger import logg
from . import thresholding

class Hierarchy:
    '''
    Class to organize a MMoCHi hierarchy. The Hierarchy is a tree with alternating subset and classification nodes for progressively annotating 
    cell types. Subset nodes define cell populations and the Hierarchy is initialized a root Subset "All", representing all events in the dataset. 
    All other Subsets originate from a Classification node. Classification nodes are defined with a list of markers (for high-confidence labeling), 
    and normally trigger selection of high-confidence events, training of a random forest classifier, and prediction. If a classification node is a 
    cutoff, it will only trigger a selection of high-confidence events, and those events will be carried forward. Subset nodes will also contain
    the cell type definitions used for high-confidence thresholding. 
    
    Initializing the Hierarchy, you can also define many classification defaults, which can be additionally customized for each Classification node. 
    
    Parameters
    ---------
    default_min_events
        The default minimum number of (or proportion of total) high-confidence events that must be identified for in order to train a random forest
        classifier with each Subset. If not enough events are identified, that Subset will be skipped.
    default_class_weight
        The default class_weight strategy for handling scoring. This is passed to sklearn.ensemble.RandomForestClassifier.
    default_clf_kwargs
        The default keyword arguments for classification. For more information about other kwargs that can be set, please see: 
        sklearn.ensemble.RandomForestClassifier In the case of batch-integrated classification, n_estimators refers to the (approximate) total trees 
        in the forest.    
    default_in_danger_noise_checker
        The default for whether to check for (and amplify or remove, respectively) in danger and noise events. In danger events are high-confidence 
        events at classification boundaries. Events labeled noise are high-confidence events whose nearest neighbors do not share the same label, 
        and are thus likely mislabeled. Can be a boolean, or "in danger only"/"noise only".
    default_is_cutoff
        The default for whether Classification nodes should be treated as a cutoff (triggering only high-confidence thresholding).
    default_features_limit
        Listlike of str or dictionary in the format {'modality_1':['gene_1','gene_2',...], 'modality_2':'All'}
        Specifies the default features allowed for training the classifier.
    default_max_training
        Specifies the default maximum number of events used for training. This directly affects training speed.
    default_force_spike_ins
        The default list of Subsets for which training events should be sampled with spike-ins from across batches, even if individual batches
        have enough events for training. This can be useful for cell types that are very heterogenous across batches.
    default_calibrate
        Default for whether to perform calibration on the prediction probabilities of the random forest classifier. Uncalibrated values reflect 
        the % of trees in agreement. Calibrated values more-closely reflect the % of calls correctly made at any given confidence level.
    load
        Either None (to initiate a new hierarchy) or a path to a hierarchy to load. Note that loading a hierarchy overrides all other defaults.
    '''
    def __init__(self, default_min_events: Union[int,float]=0.001,
                       default_class_weight: str='balanced_subsample', 
                       default_clf_kwargs: dict=dict(max_depth=20, n_estimators=100,
                                                      n_jobs=-1,bootstrap=True,verbose=True),
                       default_in_danger_noise_checker: Union[str,bool]=True, 
                       default_is_cutoff: bool=False, 
                       default_features_limit: Optional[Union[List[str],Dict[str,List[str]]]]=None,
                       default_max_training: int=20000,
                       default_force_spike_ins: List[str]=[], 
                       default_calibrate: bool=True,
                       load: Optional[str]=None):

        if not load is None:
            logg.info(f'Loading classifier from {load}...')
            self._load(load)
        else:
            self.default_min_events = default_min_events
            self.default_class_weight = default_class_weight
            self.default_in_danger_noise_checker = default_in_danger_noise_checker
            self.default_clf_kwargs = default_clf_kwargs
            self.default_features_limit = default_features_limit
            self.default_is_cutoff = default_is_cutoff
            self.default_max_training = default_max_training
            self.default_force_spike_ins = default_force_spike_ins
            self.default_calibrate = default_calibrate
            self.tree = treelib.Tree()
            self.tree.create_node("All","All",data=Subset("root")) #creating a root node
            self.reset_thresholds()
        return

    def _load(self, name: str):
        """
        Load Hierarchy as a .hierarchy
        
        Parameters
        ----------
        name
            Name of the file to load (without the .hierarchy; that will be automatically appended)
        """
        with open(name+'.hierarchy','rb') as f:
            self.__dict__ = pickle.load(f)
        logg.print('Loaded '+name+'.hierarchy')
        return
    
    def save(self, name: str):
        """
        Save Hierarchy as a .hierarchy
        
        Parameters
        ----------
        name
            Name of the file to save (without the .hierarchy; that will be automatically appended)
        """
        with open(name+'.hierarchy','wb') as f:
            pickle.dump(self.__dict__,f)
        return

    def copy(self):
        '''
        Performs a hard copy of the hierarchy (completely unlinked to the original).
        Usage: h2 = h.copy()          
            
        '''
        copied_self = Hierarchy()
        copied_self.__dict__ = pickle.loads(pickle.dumps(self.__dict__))
        return copied_self
    
    def add_classification(self, name: str, parent_subset: str,
                           markers: Union[str,List[str]], is_cutoff: Optional[bool]=None, **kwargs):
        '''
        Add a Classification beneath a Subset. Checks if it is correctly positioned in the tree before adding.
        
        Parameters
        ----------
        name
            Name of the Classification node. Names must be unique and should not contain an asterisk.
        parent_subset
            Name of the Subset that should parent this new node, use 'All' for the root node.
        markers
            The features that will be used for high-confidence thresholding to define subsets beneath this classification. During thresholding, 
            matching or similar feature names are looked up first in the provided data_key, then in the .var. See mmc.utils.marker for details 
            on how lookup can be specified.
        is_cutoff
            Whether to be a cutoff or a classifier. If None, reference the default.
        **kwargs
            arguments to pass to Classification class constructor
        '''
        assert parent_subset in self.tree.expand_tree(), f"{parent_subset} is not in hierarchy"
        assert type(self.tree[parent_subset].data) is Subset, f"{parent_subset} is not a valid subset"
        assert '*' not in name, f'name ({name}) cannot contain the following character: *'
        if not isinstance(markers, list):
            markers = [markers]
        
        # *'s placed around the id reference that this is a classification node
        self.tree.create_node("*"+name+"*",name,parent=parent_subset,data=Classification(markers,is_cutoff=is_cutoff,**kwargs))
        return
        
    def add_subset(self, name: str, parent_classification: str,
                   values: Union[List[str], List[List[str]], Dict[str,str]], 
                   color: str='#FFFFFF', **kwargs):
        '''
        Add a Subset beneath a Classification. Checks if it is correctly positioned in the tree before adding.
        
        Parameters
        ----------
        name
            Name of the Subset node. Names must be unique and should not contain an asterisk.
        parent_classification
            Name of the Classification node that should parent this new node.
        values
            list of str (e.g. ['pos','any','neg']) or list of list of strings (e.g. [['pos','any','neg'],['any','pos','neg']])
            The values used to define subsets. Lists of str should be the same length (and order) as the markers of the parent classification node.
            Lists of lists are joined as logical or's. 
            This can also be provided as a dict, which is sent to hc_defs(), allowing for more complex definitions. See hc_defs for more details.
        color
            Hex code (including the #) that defines the color associated with the subset when displaying the hierarchy or making color dictionaries
        **kwargs
            arguments to pass to Subset class constructor
        '''
        assert parent_classification in self.tree.expand_tree(), f"{parent_classification} is not in hierarchy"
        assert (type(self.tree[parent_classification].data) is Classification), f"{parent_classification} is not a valid classification layer"
        assert '*' not in name, f'name ({name}) cannot contain the following character: *'

        # Allows values to be passed as a dict of args for hc_defs
        if isinstance(values,dict):
            markers = self.tree[parent_classification].data.markers
            tag = values.copy()
            values = hc_defs(markers, **values)
        else:
            tag = None
        # Properly reformat values as a list or lists
        if not isinstance(values, list):
            values = [values]
        if not isinstance(values[0],list):
            values = [values]
        
        # Do typechecking to make sure the lists are all the proper length
        markers = self.tree[parent_classification].data.markers
        for val in values:
            assert len(markers) == len(val), "high-confidence values must match up to parent's high-confidence markers: "+str(markers)
            
        if tag is None:
            tag = self._label(name,values,markers)
        else:
            tag = self._label_dict(name,tag)
        
        # Create the node
        self.tree.create_node(tag,name, parent=parent_classification, data=Subset(values,color,**kwargs))
        return

    def _label_dict(self, name: str, values: Dict[str, str]) -> str:
        """
        Create a nicely formatted label for the node, to be used with plotting. Using a dictionary for cleaner output.
        
        Parameters
        ----------
        name
            Name of the node or classification layer
        values
             Values to add to the label, using a dictionary as input.
        Returns
        -------
        The formatted label for the node
        """
        tag = name + " ("
        if 'pos' in values:
            tag = tag + "+ ".join(values['pos']) + "+ "
        if 'neg' in values.keys():
            tag = tag + "- ".join(values['neg']) + "- "

        if 'other' in values:
            if type(values['other']) is str:
                values['other'] = [values['other']]
            tag = tag + f"=={values['OTHER']} ".join(values['other']) + f"=={values['OTHER']} "

        if 'any_of' in values:
            tag = tag[:-1] + ') and ['
            if not type(values['any_of'][0]) is list:
                values['any_of'] = [values['any_of']]
            if not 'n' in values:
                values['n'] = 1
            if not type(values['n']) is list:
                values['n'] = [values['n']]
            if (len(values['n']) == 1) and (len(values['any_of'])>1):
                values['n'] = values['n'] * len(values['any_of'])

            for any_of, n in zip(values['any_of'], values['n']):
                tag = tag + f'{n} of (' + f"+ ".join(any_of) + f"+) "

                if 'any_ofs_connector' in values and values['any_ofs_connector'] == '|':
                    tag = tag + 'or '
                    n_rm = 4
                else:
                    tag = tag + 'and '
                    n_rm = 5
            tag = tag[:-n_rm] + ']'
        else:
            tag = tag[:-1] + ')'
            
        tag = re.sub('(\((\[[1-9a-zA-Z]*(\+|\-)))\)','\g<2>',tag)
        return tag
    
    def _label(self, name: str, values: Union[List[str], List[List[str]]], markers: List[str]) -> str:
        """
        Create a nicely formatted label for the node, to be used with plotting. Using a list.

        Parameters
        ----------
        name
            Name of the node or classification layer
        values
            Values from classification subset to be added to the label
        markers
            Names of parent classification markers.
        Returns
        -------
        Returns the formatted label for the node
        """
        df = pd.DataFrame(values)
        static = []
        for col in df.columns:
            if len(df[col].unique()) == 1:
                static.append(col)
        value_dict = {'pos':'+','neg':'-'}
        if (set(list(np.array(values[0])[static]))-{'any'}) != set():
            tag = name + " ("
            for i,value in enumerate(values[0]):
                if i not in static:
                    continue
                if value in value_dict.keys():
                    value = markers[i]+value_dict[value]
                elif value !='any':
                    value = markers[i]+"=="+str(value)
                else:
                    value = ''
                tag = "".join([tag,value])
                tag = tag + " " 
                if value == '': # Remove the space if it's an any
                    tag = tag[:-1]
            if len(tag) > 0 and tag[-1] == " ":
                tag = tag[:-1]
            tag = tag + ")"
            if len(static) != len(values[0]):
                tag = tag + " and [("
        else:
            tag =  name +' [('
        for j,vals in enumerate(values):
            for i,value in enumerate(vals):
                if i in static:
                    continue
                if value in value_dict.keys():
                    value = markers[i]+value_dict[value]
                elif value !='any':
                    value = markers[i]+"=="+str(value)
                else:
                    value = ''
                tag = "".join([tag,value])
                tag = tag + " "
                if value == '': # Remove the space if it's an any
                    tag = tag[:-1]
            if len(tag) > 0 and tag[-1] == " ":
                tag = tag[:-1]
            if j != len(values)-1:
                tag = tag + ") or ("
            elif len(static) != len(values[0]):
                tag = tag + ")]"
        tag = re.sub('(\((\[[1-9a-zA-Z]*(\+|\-)))\)','\g<2>',tag)
        return tag
    
    def flatten_children(self, parent_subset_to_dissolve: str):
        '''
        Flattens child nodes of the hierarchy. The direct child classification information is moved to the same layer as the parent and their 
        parent classification is dissolved. If parents and children have conflicting definitions, edge-case behavior may lead to the generation 
        of nodes with duplicate markers (which is otherwise prevented) but should still function. This should be done before setting thresholds.
        
        Parameters
        ----------
        parent_subset_to_dissolve
            The subset that is parent to the classification layer being dissolved.
        '''
        assert parent_subset_to_dissolve in self.tree.expand_tree(parent_subset_to_dissolve), parent_subset_to_dissolve+" is not in hierarchy"
        # assert (type(self.tree[parent_subset_to_dissolve].data) is Subset), parent_subset_to_dissolve+" is not a valid subset layer"

        children = self.tree.children(parent_subset_to_dissolve)
        assert len(children) == 1, parent_subset_to_dissolve+" must have only 1 child"
        siblings = self.tree.siblings(parent_subset_to_dissolve)
        grandchildren = self.tree.children(children[0].identifier)
        grandparent = self.tree.parent(parent_subset_to_dissolve)
        parent_subset_to_dissolve = self.tree[parent_subset_to_dissolve]

        grandparent.data.markers = grandparent.data.markers+children[0].data.markers

        for grandchild in grandchildren:
            grandchild.tag = f'{grandchild.identifier}{parent_subset_to_dissolve.tag.split(parent_subset_to_dissolve.identifier)[-1]}'+ ' AND' + \
                             grandchild.tag.split(grandchild.identifier)[-1]
            grandchild.data.values = utils._list_tools(parent_subset_to_dissolve.data.values,'*',grandchild.data.values)
            self.tree.move_node(grandchild.identifier, grandparent.identifier)

        for sibling in siblings:
            sibling.data.values = utils._list_tools(sibling.data.values,'*',[['any']*len(children[0].data.markers)])

        grandparent.tag = grandparent.tag[:-1]+ "/" + children[0].tag[1:]
        self.tree.remove_node(parent_subset_to_dissolve.identifier)
        return

    def get_info(self, name: str, info_type: Union[List[str], str, dict]) -> Union[str, List[str], float, bool]:
        '''
        Gets specified information for a node in the hierarchy. Works with both Classification and Subset nodes.

        Paremters
        ---------
        name
            Node to get information from, if none gets general information on the hierarchy
        info_type
            What type of information to get for the node
        
        Returns
        -------
        info
            Information on the inputted node
        
        '''
        if type(info_type) is list:
            return [self.get_info(name,i) for i in info_type]
        assert info_type in self.tree[name].data.__dict__.keys(), f'{info_type} not in keys'
        default_info = self.__dict__['default_'+info_type]
        if name is None:
            return default_info
        info = self.tree[name].data.__dict__[info_type]
        if info is None:
            return default_info
        if type(default_info) is dict and type(info) is dict:
            info = utils._default_kwargs(info.copy(), default_info.copy())
        return info

    def get_classifications(self) -> List[str]:
        '''
        Provides all a list of all classification (or cutoff) nodes in the tree.
        '''
        nodes = []
        for node in self.tree.nodes.values():
            if type(node.data) is Classification:
                nodes.append(node.identifier)
        return nodes
    
    def subsets_info(self, name: str) -> Dict[str, List[int]]:
        '''
        Provides information of the subsets beneath a classification layer and their high-confidence threshold definitions.

        Parameters
        ----------
        name
            The name of the classification level to query

        Returns
        -------
            A dict containing each subset the given classification layer as keys and the high-confidence values associated with them.
        '''
        info = {}
        for node in self.tree.children(name):
            info[node.identifier] = node.data.values
        return info            
        
    def classification_parents(self, name: str) -> Tuple[str,str]:
        '''
        Provides the names of a node's parent and grandparent. This can be useful for subsetting.
        
        Parameters
        ----------
        name
            The name of the classification level to query

        Returns
        -------
        A tuple of strings for the parent followed by the grandparent.
        '''
        parent = self.tree.parent(name)
        if parent.identifier != "All":
            grandparent = self.tree.parent(parent.identifier)
        else:
            grandparent = parent
        return parent.identifier, grandparent.identifier
    
    def classification_markers(self, name: str) -> Tuple[List[str], Dict[str, List[int]]]:
        '''
        Provides markers used in one Classification node paired with the high-confidence definitions for each of its Subset nodes. 

        Parameters
        ----------
        name
            The name of the classification node to query

        Returns
        -------
        markers
            The list of markers used
        hc_subsets
            The list of definitions for each subset
        '''
        markers = self.tree[name].data.markers
        hc_subsets = self.subsets_info(name)
        return markers,hc_subsets
    
    def set_clf(self, name: str, clf, feature_names: List[str]):
        '''
        Stores a trained classifier and a list of features used for training of a specified classification level.
        
        Parameters
        ----------
        name
            The name of the classification node to add this information to
        clf
            A sklearn machine learning classifier to store in the hierarchy and use for prediction
        feature_names
            A list of all features used for classification
        
        '''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        assert not self.get_info(name,'is_cutoff'), name+"is a cutoff layer and cannot accept a clf" 
        self.tree[name].data.classifier = clf
        self.tree[name].data.feature_names = feature_names
        return
    
    def has_clf(self, name: str) -> bool:
        '''
        Checks whether a given node has a trained classifier defined.
        
        Parameters
        ----------
        name
            The name of the classification node to query
        Returns
        -------
        A boolean for whether the node has a classifier defined.             
        '''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        if self.tree[name].data.is_cutoff:
            return False
        return not self.tree[name].data.classifier is None
    
    def get_clf(self, name: str, base: bool=False):
        '''
        Gets the classifier and feature names of a given node. If base, returns base estimator 
        
        Parameters
        ----------
        name
            The name of the classification node to query
        base
            Whether to always return the base estimator (RandomForest) or potentially return the CalibrationCV if calibration was performed
            
        Returns
        -------
        A boolean for whether the node has a classifier defined.            
        '''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        clf = self.tree[name].data.classifier
        if base:
            if type(clf) is sklearn.calibration.CalibratedClassifierCV: 
                try: # accounts for deprecation of base_estimator
                    clf2 = clf.base_estimator
                except:
                    clf2 = 'failed'
                if type(clf2) is str:
                    clf2 = clf.estimator
                clf = clf2
        return clf, self.tree[name].data.feature_names
    
    def reset_thresholds(self):
        '''
        Creates an empty thresholds DataFrame
        '''
        self.thresholds = pd.DataFrame(columns=['minimum','maximum','interactive'])
        self.thresholds.index=pd.MultiIndex.from_frame(pd.DataFrame(columns=['marker','name','batch']))
        return
    
    def set_threshold(self, marker: str, thresholds: Tuple[float], interactive: bool,
                      name: Optional[str]=None, batch: Optional[str]=None):
        '''
        Sets a threshold in the Hierarchy for one marker
        
        Parameters
        ----------
        marker
            Marker to set threshold on
        thresholds
            List of thresholds to use for that node
        interactive
            Whether to interactively set thresholds
        name
            Name of the classification layer to define on. Use None for all.
        batch
            Batch to threshold. Use None for all.
        '''
        assert pd.isna(name) or name in self.get_classifications(), f'{name} is not a valid classification layer'
        if pd.isna(name):
            valid_markers = []
            for classification in self.get_classifications():
                valid_markers.extend(self.tree[classification].data.markers)
        else:
            valid_markers = self.tree[name].data.markers
        assert marker in valid_markers, f'{marker} is not a valid marker'

        data = [min(thresholds), max(thresholds), interactive]
        self.thresholds.sort_index(inplace=True)        
        try:
            self.thresholds.loc[(marker,name,batch)] = data
        except:
            self.thresholds = pd.concat([self.thresholds,pd.DataFrame([data],
                                         index=[(marker,name,batch)],columns=self.thresholds.columns)])
        return
    
    def batchless_thresholds(self, name: str=None, batch: str=None):
        '''
        Sets thresholds, removing any that are batch-specific, and averages those across batches
        
        Parameters
        ----------
        name
            Name of the classification layer to define on. Use None for all.
        batch
            Should be set to None.
        '''
        assert pd.isna(name) or name in self.get_classifications(), f'{name} is not a valid classification layer'
        for marker in set(list(zip(*self.thresholds.index))[0]):
            data = self.thresholds.loc[(marker,name,slice(None))].mean()
            data['interactive'] = bool(data['interactive'])
            self.thresholds = pd.concat([self.thresholds,pd.DataFrame([data],
                                         index=[(marker,name,batch)],columns=self.thresholds.columns)])
        self.thresholds.sort_index(inplace=True)
        return
    
    def drop_threshold(self, marker: str, name: Optional[str]=None, batch: Optional[str]=None):
        '''
        Remove thresholds from the database. Pass slice(None) to drop all thresholds matching the other keys.
        

        Parameters
        ----------
        marker
            Marker to drop thresholds on.
        name
            Name of the classification layer to remove thresholds for. Use slice(None) for all.
        batch
            Batch to remove thresholds for. Use slice(None) for all.
        '''
        self.thresholds.sort_index(inplace=True)
        try:
            if pd.isna(name) or pd.isna(batch):
                x, y, z = zip(*self.thresholds.index.to_flat_index())
            if pd.isna(name):
                name = pd.isna(list(y))
            if pd.isna(batch):
                batch = pd.isna(list(z))
            self.thresholds = self.thresholds.loc[self.thresholds.index.difference(
                                                  self.thresholds.loc[(marker,name,batch)].index)]
        except:
            logg.warning('No thresholds dropped')
        return
    
    def get_threshold_info(self, marker: str, name: str,
                           batch: str=None, flexible_level: bool=True,
                           flexible_batch: bool=True) -> Tuple[bool, Tuple[float, float]]:
        '''
        Identifies and returns threshold information, with support for flexibly for level or batch. 
        
        Parameters
        ----------
        marker
            Marker to find thresholding information on
        name
            Name of the node in the hierarchy
        batch
            Batch of the marker to find information on
        flexible_level
            Whether to search any level, if the specified level lacks information
        flexible_batch
            Whether to search any batch, if the specified batch lacks information
            
        Returns
        -------
        interactive: bool 
            Whether thresholds have been done interactively
        (minimum, maximum): tuple of float
            Minimum and maximum thresholds for the specified marker, node, and batch
        '''
        self.thresholds.sort_index(inplace=True)
        x=[]
        try:
            x = self.thresholds.loc[(marker,name,batch)]
        except:
            pass
        if len(x) == 0 and flexible_level:
            try:
                x = self.thresholds.loc[(marker,None,batch)]
            except:
                pass
        if len(x) == 0 and flexible_batch:
            try:
                x = self.thresholds.loc[(marker,level,None)]
            except:
                pass
        if len(x) == 0 and flexible_level and flexible_batch:
            try:
                x = self.thresholds.loc[(marker,None,None)]
            except:
                pass
        assert len(x) > 0, f'No thresholds found for {marker}, {name}, {batch}'
        if type(x) == pd.DataFrame:
            x = x.iloc[0]
        return x.interactive, (x.minimum,x.maximum)
    
    def save_thresholds(self, save_path: str=None, non_destructive: bool=True):
        '''
        Saves thresholds as a .csv file, non_desctructive saving loads in the old file and appends new definitions onto it

        Parameters
        ----------
        save_path
            Filepath where one is saving the thresholding information 
        non_destructive
            Whether to append the current file with new thresholding information or whether to overwrite any existing file
        '''
        if save_path is None:
            return self.thresholds
        else:
            df = self.thresholds.copy()
            if non_destructive:
                try:
                    df2 = pd.read_csv(save_path)
                    if not ('marker' in df.columns and 'name' in df.columns and 'batch' in df.columns):
                        df2 = df2.rename({'Unnamed: 0':'marker',
                                          'Unnamed: 1':'name',
                                          'Unnamed: 2':'batch'},axis='columns')
                    df2 = df2.set_index(['marker','name','batch'])
                    for i in df.index:
                        try:
                            df2 = df2.drop(i)
                        except:
                            pass
                    df = pd.concat([df,df2])
                    logg.info('Saving non-destructively...')
                except:
                    pass
            df.to_csv(save_path)
        return
    
    def load_thresholds(self, df: Union[str, pd.DataFrame], verbose: bool=False):
        '''
        Loads in thresholds from a `.csv` file.

        Parameters
        ----------
        df
            A reference to a DataFrame with columns for marker, name, batch, maximum, minimum, and interactive
            This should be the exact number of columns. 
        verbose
            Whether to print log messages during loading.
        '''
        if type(df) == str:
            df = pd.read_csv(df)
            if not ('marker' in df.columns and 'name' in df.columns and 'batch' in df.columns):
                df = df.rename({'Unnamed: 0':'marker',
                                'Unnamed: 1':'name',
                                'Unnamed: 2':'batch'},axis='columns')
            df = df.set_index(['marker','name','batch'])
        for index, row in df.iterrows():
            try:
                self.set_threshold(index[0], (row.maximum, row.minimum), row.interactive, index[1], index[2])
            except Exception as e:
                if verbose:
                    logg.error(f'Failed to add in thresholding at {index}')
        logg.print('Loaded thresholds.')
        return    

    def run_all_thresholds(self, adata: anndata.AnnData,
                           data_key: str=utils.DATA_KEY, batch_key: str=None,
                           mode: str='fill in', interactive: bool=True,
                           plot: bool=True, limit: Optional[Union[str, List[str]]]=None, 
                           batch_marker_order: bool=False, skip: List[str]=[]):
        '''
        Run thresholding using the `thresholding.threshold()` function. First uses `mmc.get_data` to search for marker.
        If marker is not found, gives up and ask whether to label it as interactive or not.
        
        Parameters
        ----------
        adata
            Object containing marker expression in the .X, .obsm[data_key], or .obs for high confidence thresholding.
        data_key
            obsm location to probe for other modalities beyond the .X
        batch_key
            If none, don't run with batch. Otherwise, the batch on which to threshold 
        mode
            One of: ['fill in','rerun all', 'rerun specified','every level'] 
            Whether to fill in empty holes with broad thresholds, rerun all globals, rerun all thresholds that have been
            specified (and not fill in new ones), or rerun every level separately. Add "fancy" to the mode name to trigger
            fancy plots (interactive widgets).
        interactive
            Whether to run thresholding.threshold() interactively or not
        plot
            Whether to display a plot when running thresholding.threshold()
        limit
            If specified, only thresholds marker(s) included in the limit
        batch_marker_order
            Whether to order by batch first or marker first (as the outer loop). If true, marker is the outer loop.
        skip
            Markers to skip for thresholding
        '''
        all_markers, all_levels = [], []
        for classification in self.get_classifications():
            for marker in self.tree[classification].data.markers:
                if not marker in skip:
                    all_markers.append(marker)
                    all_levels.append(classification)

        markers,levels = [],[]
        if 'fill in' in mode or 'rerun all' in mode:
            markers = list(set(all_markers))
            levels = [None]*len(markers)
        elif'rerun specified' in mode:
            to_run, levels = zip(*set((index[0,1] for index, row in self.thresholds.iterrows())))
        elif 'every level' in mode:
            markers, levels = all_markers, all_levels
        else:
            assert False,  f"{mode} is not a supported mode"
        if not limit is None:
            mask = [mark in limit for mark in markers]
            markers = list(np.array(markers)[mask])
            levels = list(np.array(levels)[mask])
            
        if 'fancy' in mode:
            fancy_resolution = []

        batch_iterable = zip(*utils.batch_iterator(adata,batch_key))
        marker_iterable = zip(markers, levels)
        if batch_marker_order:
            iterable = itertools.product(batch_iterable,marker_iterable)
        else:
            iterable = itertools.product(marker_iterable,batch_iterable)
            
        for (a,b),(c,d) in iterable:
            if batch_marker_order:
                mask, batch, marker, level = a, b, c, d  
            else:
                mask, batch, marker, level = c, d, a, b  
                
            try:
                _, t = self.get_threshold_info(marker,level,batch)
                if 'fill in' in mode:
                    continue
            except:
                t=None

            if not level is None:
                print(f'On level {level}')
            try: 
                threshold = thresholding.threshold(marker,adata[mask],data_key,preset_threshold=t,
                                         include_zeroes=False, n=0, force_model=False,
                                         plot=plot, interactive=interactive, fancy='fancy' in mode, run=False,
                                         title_addition= '' if batch_key is None else batch)
                if not 'fancy' in mode:
                    print(threshold)
                    self.set_threshold(marker,threshold, False, level, batch)
                else:
                    fancy_resolution.append((marker, threshold, level, batch))
            except Exception as e:
                # logg.warning(f'Ran into error, skipping marker...', exc_info=1)
                logg.warning(f"{marker} not found in adata, skipping thresholding.")
                    
        print('Completed!')
        if 'fancy' in mode: 
            try: 
                from ipywidgets import Button, Output
                from IPython.display import display

            except ImportError:
                raise ImportError('Please install ipywdigets using pip install ipywdigets')
            button = Button(description="Click to run thresholds!")
            output = Output()
            display(button, output)

            def on_button_clicked(b):
                with output:
                    for marker,promise, level, batch in fancy_resolution:
                        self.set_threshold(marker,thresholding._fancy_resolver(promise), False, level, batch)
                    logg.print(f'------------------------\nThreshold fancy save completed\n------------------------')

            button.on_click(on_button_clicked)

        return
    
    def get_all_markers(self) -> Set[str]:
        '''
        Provides a list of all the markers used for high-confidence thresholding.
        '''
        all_markers = []
        for name in self.get_classifications():
            all_markers.extend(self.classification_markers(name)[0])
        return set(all_markers)

    def check_all_markers(self, adata: anndata.AnnData, data_key: Optional[str]=utils.DATA_KEY):
        '''
        Asserts all markers identified by `.get_all_markers()` are in adata.X or .obsm[data_key].
        
        Parameters
        ----------
        adata
            Object containing marker expression in the .X, .obsm[data_key], or .obs for high confidence thresholding.
        data_key
            obsm location to probe for other modalities beyond the .X
        '''
        all_markers = self.get_all_markers()
        cannot_find = []
        for i in all_markers:
            try:
                utils.get_data(adata,i,data_key)
            except:
                cannot_find.append(i)
        assert not cannot_find, f'Cannot find any of {cannot_find} in adata.X or adata.obsm["{data_key}"]'
        return
    
    def color_dict(self, new_color_palette: bool = False, 
                   mode: str='LEAF_DEPTH', default_color: str='#d3d3d3',
                   **kwargs) -> Dict[str, str]:
        '''
        Provides a dictionary of colors associated with each subset in the hierarchy
        
        Parameters
        ----------
        new_color_palette
           False to return the current color_dict
           True replaces the hierarchy color system with a new cubehelix palette
        mode
           One of ["ZIGZAG", "WIDTH", "DEPTH"]
           The order to label the tree
        default_color
           Hex code of the color default, for any subsets not already defined in the color palette.
        **kwargs
           Sent to sns.cubehelix_palette, a hue of > 1, and a rot >= 1 is recommended
           
        Returns
        -------
        colors: Dict of subset of hierarchy to colors
            A dictionary of the colors for each subset of the hierarchy

        '''
        
        if self.tree.nodes:
            if 'ZIGZAG' in mode:
                all_nodes = list(self.tree.expand_tree(mode=self.tree.ZIGZAG))
            elif 'WIDTH' in mode:
                all_nodes = list(self.tree.expand_tree(mode=self.tree.WIDTH))
            elif 'DEPTH' in mode:
                all_nodes = list(self.tree.expand_tree(mode=self.tree.DEPTH))
            else:
                assert False, f'mode must be one of "ZIGZAG", "WIDTH", "DEPTH". Not: {mode}'
            if 'LEAF' in mode:
                all_nodes = [i for i in all_nodes if self.tree[i].is_leaf()]
        else:
            assert False, f'No nodes in tree'
            
        if type(new_color_palette) is list:
            assert len(new_color_palette) == len(all_nodes)
            c = new_color_palette
            new_color_palette = True
        elif type(new_color_palette) is dict:
            c = list()
            for n in all_nodes:
                c.append(new_color_palette[n] if n in new_color_palette else default_color)
            new_color_palette = True
        elif type(new_color_palette) is bool and new_color_palette:
            import seaborn as sns
            x = sns.cubehelix_palette(n_colors=len(all_nodes),**kwargs)
            c = ['#%02x%02x%02x' % tuple(int(ii*255) for ii in i) for i in x]
        
        colors = {}
        for i, n in enumerate(all_nodes):
            try:
                if new_color_palette:
                    self.tree[n].data.color = c[i]
                colors[self.tree[n].identifier] = self.tree[n].data.color
            except:
                pass
        if new_color_palette:
            return
        return colors
    
    def display(self, plot: bool=False,
                return_graph: bool=False, supress_labels: bool=False,
                node_width: int=4, node_height: int=1, 
                font_mult: float=1):
        """
        Display the tree in a user-friendly format.
        
        Parameters
        ----------
        plot
            Whether to display as a plot (True) or as text (False) 
        return_graph
            Whether to return the graph object (when plot = True)
        supress_labels
            Whether to not display the high-confidence thresholding rules used to identify subsets
        node_width
            Width of the nodes within the hierarchy display
        node_height
            Height of the nodes within the hierarchy display
        font_mult
            Font size is 16 times this font multiplier
        
        Returns
        -------
        If return_graph, returns pydot.graph_from_dot_data graph
        """
        if plot:
            import pydot
            import IPython
            graph = pydot.graph_from_dot_data(self.to_graphviz(supress_labels=supress_labels, 
                                                               node_width=node_width,
                                                               node_height=node_height,
                                                               font_mult=font_mult))[0]
            plt = IPython.display.Image(graph.create_png())
            IPython.display.display(plt)
        else:
            self.tree.show()
        if return_graph and plot:
            return graph
        return
            
    def to_graphviz(self, supress_labels: bool=False,
                    node_width: int=4, node_height: int=1,
                    font_mult: float=1) -> str:
        """
        Exports the tree in the dot format of the graphviz software, useful for plotting.
        
        Parameters
        ----------
        supress_labels
            Whether to not include high-confidence information for subsets of the hierarchy
        node_width
            Width of the nodes within the hierarchy display
        node_height
            Height of the nodes within the hierarchy display
        font_mult
            Font size is 16 times this font multiplier
            
        Returns
        -------
        string_graphviz
            Hierarchy in graphviz format for plotting 
        """
        import textwrap

        nodes, connections = [], []
        if self.tree.nodes:
            for n in self.tree.expand_tree(mode=self.tree.WIDTH):
                if n == 'All':
                    continue
                nid = self.tree[n].identifier
                if supress_labels:
                    label = nid
                else:
                    label = self.tree[n].tag
                fixedsize = True
                style='filled'
                if type(self.tree[n].data) is Classification:
                    shape = "box"
                    height = node_height
                    width = node_width
                    fontsize = 24*font_mult
                    color = '#FFFFFF'
                else:
                    if not supress_labels:
                        label = "\n(".join(label.split("(",1))
                        label = " (".join(label.split("(",1))
                        label = ") ".join(label.split(")",1))
                        label = " [".join(label.split("[",1))
                        label = "] ".join(label.split("]",1))
                        label = "\n".join(textwrap.wrap(label, width=30,break_long_words=True,replace_whitespace =False))
                        label = "(".join(label.split(" (",1))
                        label = ")".join(label.split(") ",1))
                        label = "[".join(label.split(" [",1))
                        label = "]".join(label.split("] ",1))
                        label = ")".join(label.split(") ",1))
                        label = "[".join(label.split(" [",1))
                        label = "]".join(label.split("] ",1))
                    shape = "oval"
                    height = node_height
                    width = node_width
                    fontsize = 16*font_mult
                    color = self.tree[n].data.color+'80'
                state = (f'"{nid}" [label="{label}", shape={shape}, fixedsize={fixedsize}, fontsize={fontsize},' + 
                        f'width={width}, height={height}, style={style}, fillcolor="{color}", color="black"]')
                nodes.append(state)
                for c in self.tree.children(nid):
                    cid = c.identifier
                    connections.append(f'"{nid}" -> "{cid}"')
        string_graphviz = 'digraph' + ' tree {\n'
        for n in nodes:
            string_graphviz = string_graphviz + '\t' + n + '\n'
        if len(connections) > 0:
             string_graphviz = string_graphviz + '\n'
        for c in connections:
            string_graphviz = string_graphviz + '\t' + c + '\n'
        string_graphviz = string_graphviz + '}'
        return string_graphviz   

class Subset:
    '''
    A Hierarchy building block, describing a population of cells. These can be added to a Hierarchy using the `.add_subset()` method

    Parameters
    ----------
    values
        list of str (e.g. ['pos','any','neg']) or list of list of strings (e.g. [['pos','any','neg'],['any','pos','neg']])
        The values used to define subsets. Lists of str should be the same length (and order) as the markers of the parent classification node.
        Lists of lists are joined as logical or's. 
        This can also be provided as a dict, which is sent to hc_defs(), allowing for more complex definitions. See hc_defs for more details.
    color
        Hex code (including the #) that defines the color associated with the subset when displaying the hierarchy or making color dictionaries
    '''
    def __init__(self,values: Union[List[str], List[List[str]]],color: str='#FFFFFF'):
        self.values = values
        self.color = color
        return

class Classification:
    '''
    A Hierarchy building block, describing subsetting rules. These can be added to a Hierarchy using the `.add_classification()` method

    Parameters
    ---------
    markers
        The features that will be used for high-confidence thresholding to define subsets beneath this classification. During thresholding, 
        matching or similar feature names are looked up first in the provided data_key, then in the .var. See mmc.utils.marker for details 
        on how lookup can be specified.
    min_events
        The minimum number of (or proportion of total) high-confidence events that must be identified for in order to train a random forest
        classifier with each Subset. If not enough events are identified, that Subset will be skipped.
    class_weight
        The class_weight strategy for handling scoring. This is passed to sklearn.ensemble.RandomForestClassifier. 
    in_danger_noise_checker
        Whether to check for (and amplify or remove, respectively) in danger and noise events. In danger events are high-confidence 
        events at classification boundaries. Events labeled noise are high-confidence events whose nearest neighbors do not share the same label, 
        and are thus likely mislabeled. Can be a boolean, or "in danger only"/"noise only".
    classifier
        The classifier to be used for classification. If defined, must also define feature_names.
    features_limit
        listlike of str or dictionary in the format {'modality_1':['gene_1','gene_2',...], 'modality_2':'All'}
        Specifies the default features allowed for training the classifier.
    feature_names
        Names of features used to train this classifier. Not set if classifier is None.
    is_cutoff
        The default for whether Classification nodes should be treated as a cutoff (triggering only high-confidence thresholding).
    features_limit
        Listlike of str or dictionary in the format {'modality_1':['gene_1','gene_2',...], 'modality_2':'All'}
        Specifies the default features allowed for training the classifier.
    max_training
        Specifies the default maximum number of events used for training. This directly affects training speed.
    force_spike_ins
        The default list of Subsets for which training events should be sampled with spike-ins from across batches, even if individual batches
        have enough events for training. This can be useful for cell types that are very heterogenous across batches.
    calibrate
        Default for whether to perform calibration on the prediction probabilities of the random forest classifier. Uncalibrated values reflect 
        the % of trees in agreement. Calibrated values more-closely reflect the % of calls correctly made at any given confidence level.
    clf_kwargs
        The keyword arguments for classification. For more information about other kwargs that can be set, please see: 
        sklearn.ensemble.RandomForestClassifier In the case of batch-integrated classification, n_estimators refers to the (approximate) total trees 
        in the forest.
    '''
    def __init__(self, markers: List[str],
                 min_events: Optional[Union[int, float]]=None, 
                 class_weight: Optional[Union[dict, List[dict]]]=None,
                 in_danger_noise_checker: Optional[Union[str,bool]]=None, classifier=None,
                 features_limit: Optional[List[str]]=None, 
                 feature_names: Optional[List[str]]=None,
                 is_cutoff: Optional[bool]=False, max_training: Optional[int]=None,
                 force_spike_ins=[], calibrate: Optional[bool]=None,
                 clf_kwargs: dict={}):
        self.markers = markers
        if is_cutoff == True:
            self.is_cutoff = True
        else:
            self.is_cutoff = is_cutoff
            self.min_events = min_events
            self.class_weight = class_weight
            self.in_danger_noise_checker = in_danger_noise_checker
            if not classifier is None:
                assert not feature_names is None, 'Must define features used for classification'
                self.classifier = classifier
                self.feature_names = feature_names
            else:
                self.classifier = None
                self.feature_names = None
            self.clf_kwargs = clf_kwargs
            self.features_limit = features_limit
            self.max_training = max_training
            self.force_spike_ins = force_spike_ins
            self.calibrate = calibrate
        return

def hc_defs(marker_list: List[str], pos: Union[List[str],str]=[],
            neg: Union[List[str],str]=[], other: Union[List[str],str]=[],
            POS: str='pos', NEG: str='neg', 
            OTHER: str=None, UNDEFINED: str='any',
            any_of: Union[List[str],List[List[str]]]=[], any_ofs_connector: str='&',
            n: Union[List[int],int]=1, ANY_OPTIONS: Union[List[List[str]], List[str]]=['pos','any']):
    '''
    Helper function for defining simple or complex gating strategies for high-confidence thresholding. 
    
    Parameters 
    ----------
    marker_list
        All the markers in the group, usually set by the Classification node above. Markers can also have "_lo" or "_hi" appended
        to them to specify multiple thresholds on the same marker in the same node.
    pos
        List of markers that must be positive
    neg
        List of markers that must be negative
    other
        List of markers that must be other
    POS
        Default value describing pos (should usually not be altered)
    NEG
        Default value describing neg (should usually not be altered)
    OTHER
        String or list of strings containing the value(s) to write for each marker in "other".
    UNDEFINED
        Default value describing any undefined markers (should usually not be altered)
    any_of
        Used to define more complex gating. Can be represented as a list of markers (a single grouping) or a list of these lists (multiple pairings)
    any_ofs_connector
        In the case of any_of being a list of lists, whether to join them with an "&" for 'and' or a '|' for 'or' to create complex gating, 
        "&" for "at least 1 of [CCR7,SELL] AND at least 1 of [TCF7, MAL]" 
        '|' for "at least 1 of [CD19, CD20] OR at least 2 of [JCHAIN, CD138, CD34]"
    n
        When using "any_of", how many in the group must match the threshold (e.g. n=2 is "at least 2 of"), when any_of is a list of lists, 
        this can become a list of ints the same length
    ANY_OPTIONS
        Default values describing when the condition is statisfied and alternative (should ususally not be altered).
        The alternative string is the filler when those markers are not part of the condition. 
        e.g. for "any two positive" you would use ['pos','any'] for "only two positive" you would use ['pos','neg']
        Like other options, if any_of is a list of lists, this can also become a list of lists.
    '''
    assert any_ofs_connector == '&' or any_ofs_connector == '|', f'Invalid value {any_ofs_connector} for any_ofs_connector'
    
    pos = pos if pd.api.types.is_list_like(pos) else [pos]
    neg = neg if pd.api.types.is_list_like(neg) else [neg]
    other = other if pd.api.types.is_list_like(other) else [other]
    any_of = any_of if pd.api.types.is_list_like(any_of) else [any_of]
    
    if not (any_of != [] and isinstance(any_of[0],list)):
        any_of_list = [any_of]
    else:
        any_of_list = any_of
    if not isinstance(n,list):
        ns = [n]*len(any_of_list)
    else:
        ns = n
    assert len(ns) == len(any_of_list), "Length of n's must be equal to length of any_of_list"
    any_of_markers = [item for sublist in any_of_list for item in sublist]
    all_markers = list(itertools.chain(pos,neg,any_of_markers,other))
    undefined = [mark for mark in marker_list if not mark in all_markers]
    all_markers = list(itertools.chain(pos,neg,any_of_markers,other,undefined))
    assert (set(all_markers)-set(marker_list)) == set(), f'There are missing items in all_markers: {set(all_markers)-set(marker_list)}'
    assert len(set(all_markers)) == len(all_markers), f'There are duplicates in all_markers: {sorted(all_markers)}'
    assert len(set(marker_list)) == len(marker_list), f'There are duplicates in marker_list: {sorted(marker_list)}'
    markers = []
    constants = []
    for marker in pos:
        markers.append(marker)
        constants.append(POS)
    for marker in neg:
        markers.append(marker)
        constants.append(NEG)
    for marker in undefined:
        markers.append(marker)
        constants.append(UNDEFINED)
    for i,marker in enumerate(other):
        markers.append(marker)
        if isinstance(OTHER, list):
            constants.append(OTHER[i])
        else:
            constants.append(OTHER)
    options_list_list = []
    for any_of,n in zip(any_of_list,ns):
        options_list = []
        for marker in any_of:
            markers.append(marker)
        for options in itertools.combinations(any_of,n):
            not_chosen = [mark for mark in any_of if not mark in options]
            options_list.append(hc_defs(any_of,pos=options, neg=not_chosen, POS=ANY_OPTIONS[0], NEG=ANY_OPTIONS[1])[0])
        options_list_list.append(options_list)
    
    if any_ofs_connector == '|':
        lens = [len(l) for l in any_of_list]
        for i,options_list in enumerate(options_list_list):
            anys_before = [ANY_OPTIONS[1]] * sum(lens[0:i])
            anys_after = [ANY_OPTIONS[1]] * sum(lens[i:])
            options_list_list[i] = [anys_before + options + anys_after for options in options_list]
        options_list_list = [[item for sublist in options_list_list for item in sublist]]
        
    completed_options_list = [constants]
    while len(options_list_list) > 0:
        if options_list_list[0] == []:
            options_list_list[0] = [[]]
        completed_options_list = utils._list_tools(completed_options_list, '*', options_list_list.pop(0))
    arg_sorted = [i[0] for i in sorted(enumerate(markers), key=lambda x:marker_list.index(x[1]))]
    
    sorted_options_list = [list(np.array(option)[arg_sorted]) for option in completed_options_list] 
    return sorted_options_list