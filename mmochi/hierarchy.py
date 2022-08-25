import numpy as np
import pandas as pd
import itertools
import treelib
import pickle
import gc
import re

from . import utils
from .logger import logg
from . import thresholding

class Hierarchy:
    '''
    Class to organize a hierarchy of subsets and to run and store the classifier.
    The hierarchy is a tree with alternating layers of subset nodes and classification nodes
    Subset nodes define cell populations. This hierarchy initializes with the root subset "All", 
    which represents all events in the dataset. All subsets except this root originate from a classification
    Classification nodes contain a list of markers to be used in the classification, as well as other parameters. 
    Classification nodes normally trigger  selection of ground truth events, training of a classifier, and prediction
    If a classification node is a cutoff, it will trigger only a selection of ground truth events to be carried forward.
    Subset nodes beneath classification nodes define the specific populations to be identified by the classification. 
    
    default_min_events: the fraction of events at a classification level, where the classification is skipped.
    default_class_weight: the default weighing strategy for handling scoring. This can be customized on a per-node basis.
    load: either None (to initiate a new hierarchy) or a path to a hierarchy to load, note loading a hierarchy overrides other defaults.
    '''
    def __init__(self, default_min_events=0.001, default_class_weight='balanced_subsample',load=None):
        if not load is None:
            logg.info(f'Loading classifier from {load}...')
            self._load(load)
        else:
            self.default_min_events = default_min_events
            self.default_class_weight = default_class_weight
            self.tree = treelib.Tree()
            self.tree.create_node("All","All",data=Subset("root")) #creating a root node
            self.reset_thresholds()
        return

    def _load(self,name,directory='.'):
        """Save class as a .hierarchy
        name: str, Name of the file to save
        directory: str, path to save file to, default '.' is the local directory'"""
        with open(directory+'/'+name+'.hierarchy','rb') as f:
            self.__dict__ = pickle.load(f)
        print('Loaded '+name+'.hierarchy')
        return
    
    def save(self,name,directory='.'):
        """Save class as a .hierarchy
        name: str, Name of the file to save
        directory: str, path to save file to, default '.' is the local directory'"""
        with open(directory+'/'+name+'.hierarchy','wb') as f:
            pickle.dump(self.__dict__,f)
        return

    def copy(self):
        ''' Return a hard copy of the hierarchy unlinked to the original.
        Usage: 
        h2 = h.copy()'''
        copied_self = Hierarchy()
        copied_self.__dict__ = pickle.loads(pickle.dumps(self.__dict__))
        return copied_self
    
    def add_classification(self, name,parent_subset,markers,is_cutoff=False,**kwargs):
        '''Add a classification beneath a subset, checks if it is correctly positioned in the tree
           name: Name of the classification layer
           parent_subset: Name of the subset that should parent this new node, use 'All' for the root node
           markers: The markers that will be used to define subsets beneath this classification
           is_cutoff: bool, whether to be a cutoff or a classifier node
           **kwargs: arguments to pass to Classification()'''
        assert parent_subset in self.tree.expand_tree(), parent_subset+" is not in hierarchy"
        # Check if this is a subset node or a classification layer
        assert type(self.tree[parent_subset].data) is Subset, parent_subset+" is not a valid subset"
        if not isinstance(markers, list):
            markers = [markers]
        self.tree.create_node("*"+name+"*",name,parent=parent_subset,data=Classification(markers,is_cutoff=is_cutoff,**kwargs))
        return
        
    def add_subset(self, name,parent_classification,values,color = '#FFFFFF',**kwargs):
        '''Add a subset beneath a classifier; checks it seems correctly positioned in tree, then begins to format names etc. by the definitions
        name: Name of the classification layer
        parent_classification: Name of the classification or cutoff that should parent this new node
        values: The values that will be used to define subsets beneath this classification
                   may be a list (e.g. ['pos','any','neg']) or a list of lists (e.g. [['pos','any','neg'],['any','pos','neg']])
                   it may also be a dict specifying arguments to gt_defs()'''
        # Ensure the parent is in the hierarchy and that it's a classification layer
        assert parent_classification in self.tree.expand_tree(), parent_classification+" is not in hierarchy"
        assert (type(self.tree[parent_classification].data) is Classification), parent_classification+" is not a valid classification layer"
        
        # Allows values to be passed as a dict of args for gt_defs
        if isinstance(values,dict):
            markers = self.tree[parent_classification].data.markers
            values = gt_defs(markers, **values)
            
        # Properly reformat values as a list or lists
        if not isinstance(values, list):
            values = [values]
        if not isinstance(values[0],list):
            values = [values]
        
        # Do typechecking to make sure the lists are all the proper length
        markers = self.tree[parent_classification].data.markers
        for val in values:
            assert len(markers) == len(val), "Ground truth values must match up to parent's ground truth markers: "+str(markers)
            
        # Create the node
        self.tree.create_node(self._label(name,values, markers),name,parent=parent_classification,data=Subset(values,color,**kwargs))
        return
    
    def _label(self, name,values, markers):
        """Create a nicely formatted label for the node, to be used with plotting."""
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
        tag = re.sub('(\(([[1-9a-zA-Z]*(\+|\-)))\)','\g<2>',tag)
        return tag
    
    def get_min_events(self,name):
        '''Gets the minimum events of the classifier layer. If not defined, returns the tree default.
           name: Name of the classification layer to search for
           returns 
        '''
        min_events = self.tree[name].data.min_events
        if min_events is None:
            min_events = self.default_min_events
        return min_events
    
    def get_class_weight(self,name):
        '''Returns the class weight balance of the classifier layer. If not defined, returns the tree default'''
        class_weight = self.tree[name].data.class_weight
        if class_weight is None:
            class_weight = self.default_class_weight
        return class_weight

    def get_classifications(self):
        '''Return all classification nodes in the tree'''
        nodes = []
        for node in self.tree.nodes.values():
            if type(node.data) is Classification:
                nodes.append(node.identifier)
        return nodes
    
    def subsets_info(self,name):
        '''Returns a dict containing each subset beneath a classification layer and the ground truth values associated with them'''
        info = {}
        for node in self.tree.children(name):
            info[node.identifier] = node.data.values
        return info            
        
    def classification_parents(self,name):
        '''Returns parent and grandparent information about a classifier node which can be used for subsetting.'''
        parent = self.tree.parent(name)
        if parent.identifier != "All":
            grandparent = self.tree.parent(parent.identifier)
        else:
            grandparent = parent
        return parent.identifier, grandparent.identifier
    
    def classification_markers(self,name):
        '''Returns miscellaneous information about the classifier which can be used for finding the appropriate subset.'''
        markers = self.tree[name].data.markers
        gt_subsets = self.subsets_info(name)
        return markers,gt_subsets
    
    def set_clf(self,name,clf,feature_names):
        '''Sets the classifier and feature_names of a specified classification level (name) feature_names are expected to be list-like'''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        assert not self.is_cutoff(name), name+"is a cutoff layer and cannot accept a clf" 
        self.tree[name].data.classifier = clf
        self.tree[name].data.feature_names = feature_names
        return
    
    def is_cutoff(self,name):
        '''Checks whether a given level is a cutoff layer'''
        return self.tree[name].data.is_cutoff
    
    def has_clf(self,name):
        '''Checks whether a given level has a classifier defined'''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        return not self.tree[name].data.classifier is None
    
    def get_clf(self,name):
        '''Gets the classifier and feature names of a given level'''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        return self.tree[name].data.classifier, self.tree[name].data.feature_names
    
    def get_kwargs(self,name):
        '''Get classifier kwargs defined on a given level'''
        assert type(self.tree[name].data) is Classification, name+" is not a classification layer"
        return self.tree[name].data.clf_kwargs
    
    def reset_thresholds(self):
        '''Create an empty thresholds dataframe'''
        self.thresholds = pd.DataFrame(columns=['minimum','maximum','interactive'])
        self.thresholds.index=pd.MultiIndex.from_frame(pd.DataFrame(columns=['marker','name','batch']))
        return
    
    def set_threshold(self, marker, thresholds, interactive, name=None, batch=None):
        ''' Set the threshold '''
        assert pd.isna(name) or name in self.get_classifications(), f'{name} is not a valid classification layer'
        if pd.isna(name):
            valid_markers = []
            for classification in self.get_classifications():
                valid_markers.extend(self.tree[classification].data.markers)
        else:
            valid_markers = self.tree[name].data.markers
        assert marker in valid_markers, f'{marker} is not a valid marker'

        data = [min(thresholds), max(thresholds), interactive]
        try:
            self.thresholds.loc[(marker,name,batch)] = data
        except:
            self.thresholds = pd.concat([self.thresholds,pd.DataFrame([data],index=[(marker,name,batch)],columns=self.thresholds.columns)])
        return
    
    def generate_batchless_thresholds(self, name=None, batch=None):
        ''' Set the threshold '''
        assert pd.isna(name) or name in self.get_classifications(), f'{name} is not a valid classification layer'
        for marker in set(list(zip(*self.thresholds.index))[0]):
            data = self.thresholds.loc[(marker,name,slice(None))].mean()
            data['interactive'] = bool(data['interactive'])
            self.thresholds = pd.concat([self.thresholds,pd.DataFrame([data],index=[(marker,name,batch)],columns=self.thresholds.columns)])
        return
    
    def drop_threshold(self, marker, name=None, batch=None):
        ''' Remove thresholds from the database
        pass slice(None) to drop all thresholds matching the other keys'''
        try:
            if pd.isna(name) or pd.isna(batch):
                x, y, z = zip(*self.thresholds.index.to_flat_index())
            if pd.isna(name):
                name = pd.isna(list(y))
            if pd.isna(batch):
                batch = pd.isna(list(z))
            self.thresholds = self.thresholds.loc[self.thresholds.index.difference(self.thresholds.loc[(marker,name,batch)].index)]
        except:
            logg.warning('No thresholds dropped')
        return
    
    def get_threshold_info(self, marker, name, batch=None, flexible_level=True, flexible_batch=False):
        '''Identifies threshold information (flexibly for level or batch) TODO, write more'''
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
    
    def save_thresholds(self,save_path=None,non_destructive=True):
        '''Saves thresholds as a .csv file, non_desctructive saving loads in the old file and appends new definitions onto it'''
        if save_path is None:
            return self.thresholds
        else:
            df = self.thresholds.copy()
            if non_destructive:
                try:
                    df2 = pd.read_csv(save_path).set_index(['marker','name','batch'])
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
    
    def load_thresholds(self,df=None,verbose=False):
        '''Loads in thresholds from a .csv file, Expects columns for marker, name, batch, maximum, minimum, and interactive'''
        if type(df) == str:
            df = pd.read_csv(df).set_index(['marker','name','batch'])
        for index, row in df.iterrows():
            try:
                self.set_threshold(index[0], (row.maximum, row.minimum), row.interactive, index[1], index[2])
            except Exception as e:
                if verbose:
                    logg.error(f'Failed to add in thresholding at {index}')#,exc_info=1)
        print('Loaded thresholds.')
        return    

    def run_all_thresholds(self,adata,data_key='protein',batch_key=None,mode='fill in',interactive=True,plot=True,plot_all=True,limit=None,batch_marker_order=False):
        ''' Run thresholding using the thresholding.threshold() function. First search in the data_key, 
        then remove _gex and search the .X, then give up and ask whether to label it as interactive or not.
        adata: anndata object
        data_key: obsm location to probe for other modalities beyond the .X
        batch_key: If none, don't run with batch. Otherwise, a str
        mode: ['fill in','rerun all', 'rerun specified','every level'] whether to fill in empty holes 
              with broad thresholds, rerun all globals, rerun all thresholds that have been specified 
              (and not fill in new ones), or rerun every level separately
              Add "fancy" to the mode name to trigger fancy plots
        interactive: whether to run thresholding.threshold() interactively or not
        plot: whether to display a plot when running thresholding.threshold()
        plot_all: whether to plot those you skip when doing fill in
        limit: whether to specify a specific marker or list of marker
        batch_marker_order: whether to order by batch first or marker first (as the outer loop)
        TODO, make sure this batch thing prints every time correctly
        '''
        all_markers, all_levels = [], []
        for classification in self.get_classifications():
            for marker in self.tree[classification].data.markers:
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
                logg.warning(f'Ran into error, skipping marker...', exc_info=1)
                input(f"Marker {marker} not found in adata, press enter to continue. ")
                    
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
                    print(f'------------------------\nSave completed\n------------------------')

            button.on_click(on_button_clicked)

        return
    
    def get_all_markers(self):
        all_markers = []
        for name in self.get_classifications():
            all_markers.extend(self.classification_markers(name)[0])
        return set(all_markers)

    def check_all_markers(self, adata, data_key):
        all_markers = self.get_all_markers()
        cannot_find = []
        for i in all_markers:
            try:
                utils.get_data(adata,i,data_key)
            except:
                cannot_find.append(i)
        assert not cannot_find, f'Cannot find any of {cannot_find} in adata.X or adata.obsm["{data_key}"]'
        return
    
    def color_dict(self, new_color_palette = False, mode='ZIGZAG', **kwargs):
        '''returns a dictionary of colors associated with each subset in the hierarchy
           new_color_palette: False to return the current color_dict, True replaces the hierarchy color system with a new cubehelix palette
           mode: one of "ZIGZAG", "WIDTH", "DEPTH", for the order to label the tree
           **kwargs: sent to sns.cubehelix_palette, I recommend a hue of > 1, and a rot >= 1'''
        if self.tree.nodes:
            if mode == 'ZIGZAG':
                all_nodes = list(self.tree.expand_tree(mode=self.tree.ZIGZAG))
            elif mode == 'WIDTH':
                all_nodes = list(self.tree.expand_tree(mode=self.tree.WIDTH))
            elif mode == 'DEPTH':
                all_nodes = list(self.tree.expand_tree(mode=self.tree.DEPTH))
            else:
                assert False, f'mode must be one of "ZIGZAG", "WIDTH", "DEPTH". Not: {mode}'
            if new_color_palette:
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
    
    def display(self,plot=False):
        """Display the tree in a user-friendly format.
        
        If plot is True, requires textwrap and pydot. This can be done with:        
        pip install textwrap
        pip install pydot
        sudo apt-get install graphviz
        """
        if plot:
            import pydot
            import IPython
            graph = pydot.graph_from_dot_data(self.to_graphviz())[0]
            plt = IPython.display.Image(graph.create_png())
            IPython.display.display(plt)
        else:
            self.tree.show()
        return
            
    def to_graphviz(self):
        """Exports the tree in the dot format of the graphviz software, useful for plotting."""
        import textwrap

        nodes, connections = [], []
        if self.tree.nodes:
            for n in self.tree.expand_tree(mode=self.tree.WIDTH):
                if n == 'All':
                    continue
                nid = self.tree[n].identifier
                label = self.tree[n].tag
                fixedsize = True
                style='filled'
                if type(self.tree[n].data) is Classification:
                    shape = "oval"
                    height = 1
                    width = 4
                    fontsize = 24
                    color = '#FFFFFF'
                else:
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
                    shape = "box"
                    height = 1
                    width = 3.5
                    fontsize = 16
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
    ''' Basic building block of the hierarchy. Describes a population of cells. 
    values: list the length of the markers of the Classification or Cutoff above, 
               or list of lists of 
    color: hex code (including the #) that defines the color of the subsets to be used for making color dictonaries later
        values may be a list of lists or a list, color must be a hex code'''
    def __init__(self,values,color='#FFFFFF'):
        self.values = values
        self.color = color
        return

class Classification:
    ''' Class used to build the hierarchy. 
    is_cutoff: bool, whether to be a cutoff or a classifier node
    markers: Ground truth markers, list of features used to define ground truth subsets
    min_events: None defers to default of the hierarchy. Otherwise, an int or fraction. 
    class_weight: None defers to default of the hierarchy. Otherwise, a valid class_weight in sklearn.ensemble.RandomForestClassifier
    classifier: A stored classifier 
    feature_names: Names of features used in this classifier
    Use custom min_events here to override the fraction in the overall hierarchy'''
    def __init__(self,markers, min_events=None,class_weight=None,classifier=None, feature_names = None, is_cutoff=False,clf_kwargs={}):
        self.markers = markers
        if is_cutoff:
            self.is_cutoff = True
        else:
            self.is_cutoff = False
            self.min_events = min_events
            self.class_weight = class_weight
            self.classifier = classifier
            self.clf_kwargs = clf_kwargs
            self.feature_names = feature_names
        return

def gt_defs(marker_list, pos=[], neg=[],other=[], any_of=[],n=1,OTHER=None,POS='pos',NEG='neg',UNDEFINED='any',ANY_OPTIONS=['pos','any'],any_ofs_connector='&'):
    ''' Helper function for defining simple or complex gating strategies for ground truth. 
    
    marker_list: all the markers in the group, usually set by the classification/cutoff step above
    pos, neg, or other: lists of markers that must be positive, negative, or other
    POS, NEG, or OTHER: strings that contain the value to write for those groups, OTHER may also be passed as a list of strings the same length as other. 
    UNDEFINED: string that contains the value to write for any undefined groups
    
    any_of: Used to define more complex gating. Can be represented as a list of markers (a single grouping) or a list of these lists (multiple pairings)
    any_ofs_connector: In the case of any_of being a list of list, whether to join them with an "&" for 'and' or a '|' for 'or' to create complex gating like —
                       "&" for "at least one of [CCR7,SELL] and at least one of [TCF7, MAL]" 
                       '|' for "any one of [CD19, CD20] OR any 2 of [JCHAIN, CD138, CD34]"
    
    n: int, how many in the group must be pos (e.g. n=2 is "any 2 positive"), when any_of is a list of lists, this can become a list of ints the same length

    ANY_OPTIONS is a list of two string, the value and alternative. The value string is put in place of where the any condition is satisfied, and the
    alternative string is the filler when those markers are not part of the condition. 
    e.g. for "any two positive" you would use ['pos','any'] for "only two positive" you would use ['pos','neg']
    Like other options, if any_of is a list of lists, this can also become a list of lists.
    
    todo, add examples
    todo, give this a label making function to beautify labels I create later with the knowledge this function has now
    '''
    assert any_ofs_connector == '&' or any_ofs_connector == '|', f'Invalid value {any_ofs_connector} for any_ofs_connector'
    if not pd.api.types.is_list_like(pos):
        pos=[pos]
    if not pd.api.types.is_list_like(neg):
        neg=[neg]
    if not pd.api.types.is_list_like(other):
        other=[other]
    if not pd.api.types.is_list_like(any_of):
        any_of=[any_of]
    
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
    assert sorted(all_markers) == sorted(marker_list), f'Not all markers used are in marker_list, or there are duplicate definitions for markers within {sorted(all_markers)}'
    assert len(set(marker_list)) == len(marker_list), 'There are duplicates in marker_list'
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
            options_list.append(gt_defs(any_of,pos=options, neg=not_chosen, POS=ANY_OPTIONS[0], NEG=ANY_OPTIONS[1])[0])
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
        completed_options_list = utils.list_tools(completed_options_list, '*', options_list_list.pop(0))
    arg_sorted = [i[0] for i in sorted(enumerate(markers), key=lambda x:marker_list.index(x[1]))]
    
    sorted_options_list = [list(np.array(option)[arg_sorted]) for option in completed_options_list] 
    return sorted_options_list


