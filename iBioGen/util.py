
import collections
import glob
import functools
import logging
import numpy as np
import os
import pandas as pd
import shlex
import subprocess
import sys

from sklearn.linear_model import LinearRegression

LOGGER = logging.getLogger(__name__)

## Custom exception class
class iBioGenError(Exception):
    """ General iBioGen exception """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


def detect_cpus():
    """
    Detects the number of CPUs on a system. This is better than asking
    ipyparallel since ipp has to wait for Engines to spin up.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if "NUMBER_OF_PROCESSORS" in os.environ:
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    return 1 # Default


def progressbar(nsims, finished, msg=""):
    """ prints a progress bar """
    progress = 100*(finished / float(nsims))
    hashes = '#'*int(progress/5.)
    nohash = ' '*int(20-len(hashes))
    print("\r  [{}] {:>3}% {} ".format(hashes+nohash, int(progress), msg), end="")
    sys.stdout.flush()


def _expander(namepath):
    """ expand ./ ~ and ../ designators in location names """
    if "~" in namepath:
        namepath = os.path.expanduser(namepath)
    else:
        namepath = os.path.abspath(namepath)
    return namepath


def sample_param_range(param, nsamps=1, loguniform=False):
    """ Sample a parameter from a range. This is used to allow parameter
    values in the params file to be specified as a tuple and have simulations
    sample from the range on the fly.
    """
    LOGGER.debug("Sampled from range {}".format(param))
    if isinstance(param, tuple):
        if isinstance(param[0], float):
            if loguniform:
                param = np.random.uniform(np.log10(param[0]), np.log10(param[1]), size=nsamps)
                param = np.power(10, param)
            else:
                param = np.round(np.random.uniform(param[0], param[1], nsamps), 5)
        else:
            if loguniform:
                param = np.random.uniform(np.log10(param[0]), np.log10(param[1]), nsamps)
                param = np.int32(np.power(10, param))
            else:
                param = np.random.randint(param[0], param[1], nsamps)
    elif param == 0:
        param = np.round(np.random.random(nsamps), 3)
    else:
        param = [param] * nsamps
    LOGGER.debug("Sampled Value: {}".format(param))
    return param


def tuplecheck(newvalue, dtype=str):
    """
    Takes a string argument and returns value as a tuple.
    Needed for paramfile conversion from CLI to set_params args
    """

    ## If it's a list then this is probably api mode so the types
    ## of the values should be fine.
    if isinstance(newvalue, tuple):
        ## Already a tuple so we good to go, just make sure the
        ## values are the right dtype
        newvalue = (dtype(newvalue[0]), dtype(newvalue[1]))
    elif isinstance(newvalue, list):
        try:
            newvalue = tuple(newvalue)
        except TypeError:
            raise iBioGenError("tuplecheck failed for {}, improper list format".format(newvalue))
    else:
        try:
            ## If it's just one value of the proper dtype this should
            ## suffice to catch it.
            newvalue = dtype(newvalue)
        except Exception as inst:
            ## Failed to cast to dtype, so this is probably a prior range
            ## In some cases you may get string representations of float
            ## values that _should_ be int. Here we first cast the string
            ## to a float, and then cast to dtype.
            ## https://stackoverflow.com/questions/1841565/valueerror-invalid-literal-for-int-with-base-10
            try:
                newvalue = newvalue.rstrip(")").strip("(")
                minval = dtype(float(newvalue.split("-")[0].strip()))
                maxval = dtype(float(newvalue.split("-")[1].strip()))
                newvalue = tuple([minval, maxval])
            except Exception as inst:
                raise iBioGenError("{}\ttuplecheck() failed to cast to {} - {}"\
                            .format(inst, dtype, newvalue))

    LOGGER.debug("Returning {} - {}".format(type(newvalue), newvalue))
    return newvalue


def set_params(data, param, newvalue, quiet=True):
    """
    Set a parameter to a new value. Raises error if newvalue is wrong type.
    This is used to set parameters on both the Region and LocalCommunity
    paramsdicts.

    Parameters
    ----------
    param : str
        The string name (e.g., "project_dir") for the parameter 
        that will be changed.
    
    newvalue : int, str, or tuple
        The new value for the parameter selected for `param`.
        If the wrong type is entered for newvalue
        (e.g., a str when it should be an int), an error will be raised.
        Further information about each parameter is also available
        in the documentation.
    """
    LOGGER.debug("set param: {} {} = {}".format(data, param, newvalue))
    #allowed_params = list(data.paramsdict.keys())
    ## require parameter recognition
    if not param in list(data.paramsdict.keys()):
        raise iBioGenError("Parameter key not recognized: {}"\
                                .format(param))
    try:
        data._paramschecker(param, newvalue, quiet)
    except Exception as inst:
        raise iBioGenError(BAD_PARAMETER.format(param, inst, newvalue))
    return data


def load_sims(sims, sep=" ", nrows=None, calc_slopes=True):
    """
    Load simulations either from a dataframe or from a file.

    nrows:          Number of simulations to load. If None then loads all.
    calc_slopes:    Calculate slope of pi/speciation rate regression and
                    whether or not the slope is significantly different
                    from zero using a permutation approach.
    """

    if isinstance(sims, str):
        sim_df = pd.read_csv(sims, sep=sep, header=0, nrows=nrows)
        # Get the nicely formatted params and stats
        # Carve off the last two elements which are data and the tree
        params_df = sim_df.loc[:, sim_df.columns[:-2]]

        sims = []
        for rec in sim_df["data"]:
        # split the records for each species, separated by ','
            dat = rec.split(",")
            dat = {x:{"abundance":int(y), "pi":float(z), "r":float(aa), "lambda_":float(bb)} for x, y, z, aa, bb in map(lambda x: x.split(":"), dat)}
            sims.append(dat)
        dat_df = pd.DataFrame(sims)

        # Calculate ClaDS m (a compound parameter) as a convenience
        params_df["ClaDS_m"] = params_df["ClaDS_alpha"] *\
                                    np.exp(params_df["ClaDS_sigma"]**2/2)

        # Hill numbers for pi and abundance
        for i in range(1, 4):
            params_df[f"abund_h{i}"] = dat_df.apply(_hill_number, order=i, axis=1)
        for i in range(1, 4):
            params_df[f"pi_h{i}"] = dat_df.apply(_hill_number, vals="pi", order=i, axis=1)
        for i in range(1, 4):
            params_df[f"sp_h{i}"] = dat_df.apply(_hill_number, vals="sp", order=i, axis=1)

        if calc_slopes:
            params_df["pi_lambda_slope"] = dat_df.apply(_slope, axis=1)
            params_df["slope_sign"] = dat_df.apply(_test_significance, axis=1)
    else:
        raise iBioGenError("Input simulations not understood. Must be a file name or a pandas DataFrame")
    return params_df, dat_df, sim_df["tree"]


def _slope(data):
    """
    Calculate the slope of the regression between pi and speciation rate.
    Input is one row of data from a simulation
    """
    pis = np.array([x["pi"] for x in data])
    lambdas = np.array([x["lambda_"] for x in data])
    rgr = LinearRegression().fit(pis.reshape(-1, 1),
                                 lambdas.reshape(-1, 1))
    return rgr.coef_[0][0]


def _test_significance(data, replicates=100):
    rng = np.random.default_rng()
    sl = _slope(data)

    pis = np.array([x["pi"] for x in data])
    lambdas = np.array([x["lambda_"] for x in data])

    reps = []
    for rep in range(replicates):
        rng.shuffle(pis)
        rgr = LinearRegression().fit(pis.reshape(-1, 1),
                                 lambdas.reshape(-1, 1))
        reps.append(rgr.coef_[0][0])
    lower, upper = np.quantile(reps, [0.025, 0.975])

    if sl < lower:
        return "negative"
    if sl > upper:
        return "positive"
    else:
        return "zero"


def _hill_number(data, vals="", order=1):
    """
    Calculate hill numbers from a row of simulation data

    :param array-like data: A row of simulated data from the output file
    :param vals str: Must be one of "pi", "sp", or "". The values to calculate
                        hill numbers for. If blank then just abundance.
    :param int order: The order of the hill number to calculate
    """

    if vals == "pi":
        vals = np.array([x["pi"] for x in data])
    elif vals == "sp":
        vals = np.array([x["lambda_"] for x in data])
    else:
        vals = None

    abundance = np.array([x["abundance"] for x in data])

    return _generalized_hill_number(abundance, vals=vals, order=order)


def _generalized_hill_number(abunds, vals=None, order=1, scale=True, verbose=False):
    """
    This is the Chao et al (2014) generalized Hill # formula. Generalized
    function to calculate one Hill number from a distribution of values of
    some statistic and abundances.

    :param array-like abunds: An `array-like` of abundances per species.
    :param array-like vals: An `array-like` of values per species (e.g.
        pi values per species, or trait values per species). If this parameter
        is empty then only abundance Hill numbers are calculated.
    :param float order: The Hill number to calculate. 0 is species richness.
        Positive values calculate Hill numbers placing increasing weight on
        the most abundant species. Negative values can also be specified
        (placing more weight on the rare species), but these are uncommonly
        used in practice.
    :param bool scale: Whether to scale to effective numbers of species, or
        return the raw attribute diversity. Equivalent to equation 5c in
        Chao et al 2014. You will almost never want to disable this.

    :return float: The generalized Hill number of order `order` for the given
        data axis using the formula proposed by Chao et al (2014).
    """
    ## Degenerate edge cases can cause all zero values, particulary for pi
    ## in which case we bail out immediately
    if not np.any(abunds):
        return 0

    ## Be sure abundance is scaled to relative abundance and convert to np
    abunds = np.array(abunds)/np.sum(abunds)

    ## If vals is empty then populate the vector with a list of ones
    ## and this function collapses to the standard Hill number for abundance
    if vals is None:
        vals = np.ones(len(abunds))

    ## Make sure vals is an np array or else order > 2 will act crazy
    vals = np.array(vals)
    if verbose: print(("sums:", "dij", np.sum(vals), "pij", np.sum(abunds)))
    ## sum of values weighted by abundance
    V_bar = np.sum(vals*abunds)
    if verbose: print(("vbar", V_bar))

    ## Use the special formula for order = 1
    if order == 1:
        proportions = vals*(abunds/V_bar)
        h = np.exp(-np.sum(proportions * np.log(abunds/V_bar)))
    else:
        h = np.sum(vals*(abunds/V_bar)**order)**(1./(1-order))
    if scale: h = h/V_bar
    return h


# all branches of given tree will be rescaled to TARGET_AVG_BL
TARGET_AVG_BL = 1

def _add_dist_to_root(tre):
    """
    Add distance to root (dist_to_root) attribute to each node
    :param tre: ete3.Tree, tree on which the dist_to_root should be added
    :return: void, modifies the original tree
    """

    for node in tre.traverse("preorder"):
        if node.is_root():
            node.add_feature("dist_to_root", 0)
        elif node.is_leaf():
            node.add_feature("dist_to_root", getattr(node.up, "dist_to_root") + node.dist)
            # tips_dist.append(getattr(node.up, "dist_to_root") + node.dist)
        else:
            node.add_feature("dist_to_root", getattr(node.up, "dist_to_root") + node.dist)
            # int_nodes_dist.append(getattr(node.up, "dist_to_root") + node.dist)
    return None


def _rescale_tree(tre, target_avg_length):
    """
    Returns branch length metrics (all branches taken into account and external only)
    :param tre: ete3.Tree, tree on which these metrics are computed
    :param target_avg_length: float, the average branch length to which we want to rescale the tree
    :return: float, resc_factor
    """
    # branch lengths
    dist_all = [node.dist for node in tre.traverse("levelorder")]

    all_bl_mean = np.mean(dist_all)

    resc_factor = all_bl_mean/target_avg_length

    for node in tre.traverse():
        node.dist = node.dist/resc_factor

    return resc_factor


def CBLV(tree_input, sampling_proba):
    """Rescales all trees from tree_file so that mean branch length is 1,
    then encodes them into full tree representation (most recent version)
    Compact Bijective Ladderized Vector (CBLV) from Voznika et al 2021

    :param tree_input: ete3.Tree, that we will represent in the form of a vector
    :param sampling_proba: float, value between 0 and 1, presumed sampling probability value
    :return: pd.Dataframe, encoded rescaled input trees in the form of most recent, last column being
     the rescale factor
    """

    def real_polytomies(tre):
        """
        Replaces internal nodes of zero length with real polytomies.
        :param tre: ete3.Tree, the tree to be modified
        :return: void, modifies the original tree
        """
        for nod in tre.traverse("postorder"):
            if not nod.is_leaf() and not nod.is_root():
                if nod.dist == 0:
                    for child in nod.children:
                        nod.up.add_child(child)
                    nod.up.remove_child(nod)
        return

    def get_not_visited_anc(leaf):
        while getattr(leaf, "visited", 0) >= len(leaf.children)-1:
            leaf = leaf.up
            if leaf is None:
                break
        return leaf

    def get_deepest_not_visited_tip(anc):
        max_dist = -1
        tip = None
        for leaf in anc:
            if leaf.visited == 0:
                tip = leaf
        return tip

    def get_dist_to_root(anc):
        dist_to_root = getattr(anc, "dist_to_root")
        return dist_to_root

    def get_dist_to_anc(feuille, anc):
        dist_to_anc = getattr(feuille, "dist_to_root") - getattr(anc, "dist_to_root")
        return dist_to_anc

    def encode(anc):
        leaf = get_deepest_not_visited_tip(anc)
        #print("{}".format(leaf.name), end="\t")
        yield get_dist_to_anc(leaf, anc)
        leaf.visited += 1
        anc = get_not_visited_anc(leaf)

        if anc is None:
            return
        anc.visited += 1
        yield get_dist_to_root(anc)
        for _ in encode(anc):
            yield _

    def complete_coding(encoding, max_length):
        add_vect = np.repeat(0, max_length - len(encoding))
        add_vect = list(add_vect)
        encoding.extend(add_vect)
        return encoding

    def refactor_to_final_shape(result_v, sampling_p, max_length):
        tips_coor = np.arange(0, max_length, 2)
        int_nodes_coor = np.arange(1, max_length + 1, 2)

        reshape_coordinates = np.append(int_nodes_coor, tips_coor)
    
        # reorder the columns
        result_v = result_v.iloc[:,reshape_coordinates]

        return result_v

    # local copy of input tree
    tree = tree_input.copy().treenode

    if len(tree) <= 200:
        max_len = 400
    else:
        max_len = 999

    # remove the edge above root if there is one
    if len(tree.children) < 2:
        tree = tree.children[0]
        tree.detach()

    # set to real polytomy
    real_polytomies(tree)

    # rescale branch lengths
    rescale_factor = _rescale_tree(tree, target_avg_length=TARGET_AVG_BL)

    # set all nodes to non visited:
    for node in tree.traverse():
        setattr(node, "visited", 0)

    _add_dist_to_root(tree)

    tree_embedding = list(encode(tree))
    tree_embedding = complete_coding(tree_embedding, max_len)
    #tree_embedding.append(rescale_factor)

    result = pd.DataFrame(tree_embedding, columns=[0])

    result = result.T
    # refactor to final shape: add sampling probability, put features in order

    result = refactor_to_final_shape(result, sampling_proba, max_len)

    return tree, result, rescale_factor


## Error messages
BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """

if __name__ == "__main__":
    print("Watdo")
