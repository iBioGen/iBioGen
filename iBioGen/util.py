
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


def load_sims(sims, sep=" ", calc_slopes=True):
    """
    Load simulations either from a dataframe or from a file.

    calc_slopes:    Calculate slope of pi/speciation rate regression and
                    whether or not the slope is significantly different
                    from zero using a permutation approach.
    """
    if isinstance(sims, str):
        sim_df = pd.read_csv(sims, sep=sep, header=0)
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
    return params_df, dat_df


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


## Error messages
BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """

if __name__ == "__main__":
    print("Watdo")
