
# iBioGen
Software for the study of island biodiversity dynamics from genomic data.

The development of this package was supported by the iBioGen project, funded
by the European Union’s Horizon 2020 research and innovation programme under
grant agreement No 810729. For more information please visit the project
website: https://www.ibiogen.eu/

## Major Components
The iBioGen package is composed of a suite of tools for studying ecological,
evolutionary, and multi-scale processes within island communities. Major components
of the package include:

### A joint model of macroevolution and abundance/genetic diversity
A model for studying diversification patterns, abundance distributions, and
community-scale genetic data of island radiations. The model predicts phylogenies
under a birth/death process with abundances evolved along the branches, and
genetic diversity calculated based on historical changes of abundance through time.
Abundance can evolve in one of two ways: 1) as a BM process where abundance (n) after
time (dt) is sampled from a normal distribution with mean n and standard deviation 
sigma * dt (where sigma is a model parameter), or 2) the rate of change (r) of abundance
can evolve as BM in which case abundance (n) changes through time (dt) via n * exp(r\*dt).
The speciation rate can also shift at branching events in the manner of ClaDS (Maliet et al 2019).
In either case abundances are split by random fission at speciation events (ancestral abundance
is randomly split between the two descendant lineages)

This is a stochastic simulation model, so parameters (or parameter ranges) are specified
and numerous simulations may be generated. For each simulation the model outputs a dated
phylogenetic tree, and abundance, genetic diversity, and trait values for each species
(currently only neutral traits are supported). In the simplest case simulations can be
used to characterize joint patterns in the data under various parameterizations.
Simulations can also be used to estimate model parameters from empirical data, for
example by using approximate Bayesian computation or a machine learning approach.

## Distribution
The iBioGen package is distributed in a couple different formats, so please
choose what serves your needs best:
 * [R package](#r-package)
 * [CLI/Python API](#command-line-standalone-version)

## R Package

### R Package Installation
The iBioGen R package can be installed using `devtools::install_github` while
we navigate the CRAN submission process.

```
install.packages("devtools")
library(devtools)
install_github("iBioGen/iBioGen")
library(iBioGen)
```

### R Package Usage
The iBioGen R package currently has one very feature rich function `sim_tree()`.
This function takes numerous arguments to parameterize the simulations, and
returns a matrix of simulation results, with simulations on the rows, and parameters
and data generated on the columns (for more details see the **[Output](#output)**
section below). Here we demonstrate a few different simulation parameters (you may
invoke `?sim_tree` in the normal way for detailed help info):

```
# Generate 5 simulations and store the results in `res`
res = iBioGen::sim_tree(nsims=5)

# Generate 1 simulation with 200 tips
res = iBioGen::sim_tree(ntaxa=200)

# Simulate 1 tree with 300 tips using birth_rate of 2 with the ClaDS model enabled
res = iBioGen::sim_tree(birth_rate=2, ClaDS=TRUE, ClaDS_alpha=0.7, ntaxa=300)
```
**NB:** The first time you run `sim_tree()` there is a little housekeeping which
can take ~5 minutes. This only happens once.

## Command line standalone version
The iBioGen package also exists as a standalone command line utility written
in python. The CLI provides a number of benefits including massive parallelization.

### Installation

* Install [conda](https://docs.conda.io/en/latest/miniconda.html) for python3
* `conda create -n iBioGen python=3.7`
* `conda activate iBioGen`
* `conda install -c conda-forge -c iovercast ibiogen`

### Command Line Usage
Create a params file:

    iBioGen -n wat

Look at the params and edit them if you wish:

    ------- iBioGen params file (v.0.0.8)-------------------------------------------
    watdo                ## [0] [simulation_name]: The name of this simulation scenario
    ./default_iBioGen    ## [1] [project_dir]: Where to save files
    1                    ## [2] [birth_rate]: Speciation rate
    taxa                 ## [3] [stop_criterion]: Whether to stop on ntaxa or time
    20                   ## [4] [ntaxa]: Number of taxa to simulate if stop is `ntaxa`
    4                    ## [5] [time]: Amount of time to simulate if stop is `time`
    abundance            ## [6] [process]: Whether to evolve `abundance` or growth `rate` via BM
    True                 ## [7] [ClaDS]: Whether to allow speciation rates to change along the branches a la ClaDS
    50000                ## [8] [abundance_mean]: Ancestral abundance at time 0
    0.1                  ## [9] [abundance_sigma]: Rate at which abundance changes if process is `abundance`
    0                    ## [10] [growth_rate_mean]: Ancestral population growth rate at time 0.
    0.01                 ## [11] [growth_rate_sigma]: Rate at which growth rate changes if process is `rate`
    0.1                  ## [12] [ClaDS_sigma]: Rate at which speciation rate changes if ClaDS is True
    0.9                  ## [13] [ClaDS_alpha]: Rate shift if ClaDS is True
    500                  ## [14] [sequence_length]: Length of the genomic region simulated, in base pairs
    1e-05                ## [15] [mutation_rate]: Mutation rate per base per generation
    10                   ## [16] [sample_size]: Number of samples to draw for calculating genetic diversity
    None                 ## [17] [abundance_scaling]: Scaling abundance to Ne. Can be None, log, ln or a ratio

Run 10 simulations:

    iBioGen -p params-wat.txt -s 10

Run 10 simulations on 10 cores in parallel:

    iBioGen -p params-wat.txt -s 10 -c 10

### Output
Results are written to `<project_dir>/<simulation_name>` (so in the example:
`default_iBioGen/wat-SIMOUT.csv`). Not generally human readable the results
file contains the parameters used to generate each simulation, as well as the
observed numbers of tips, the observed simulation time, and the calculated
extinction rate (as a fraction of birth events), data for each tip, including
abundance, genetic diversity, growth rate, and speciation rate, and finally
the dated tree in newick form. Field names are as follows (along with 1
example simulation):

    birth_rate stop_criterion ntaxa time process ClaDS abundance_mean abundance_sigma growth_rate_mean growth_rate_sigma ClaDS_sigma ClaDS_alpha sequence_length mutation_rate sample_size abundance_scaling obs_ntaxa obs_time turnover_rate data tree
    1 taxa 20 4 abundance False 50000 0.1 0 0.01 0.1 0.1 500 1e-05 10 None 20 3.517589203795064 0.09523809523809523 r19:5277:0.20435555555555526:0.011195997657183444:1,r18:2348:0.05888888888888879:0.010817643664283555:1,r17:3445:0.1972000000000004:0.00976585185590869:1,r16:977:0.11759999999999997:0.00956274991175375:1,r15:61:0.0031999999999999997:0.011362725225983333:1,r14:16539:0.4827111111111095:0.017471769225183797:1,r13:23182:0.6710666666666577:0.021756148739551756:1,r12:348:0.005466666666666666:0.02736456054653283:1,r11:515:0.0124:0.03659168443387108:1,r10:2466:0.04511111111111108:0.03861433944456937:1,r9:52:0.001911111111111111:0.04392777318481233:1,r8:357:0.009111111111111113:0.027613860772168673:1,r7:15:0.0:0.025752062822845163:1,r6:36:0.0008:0.02570337859771903:1,r5:4:0.0:0.02775490802853657:1,r4:402:0.01315555555555556:0.03211086569690778:1,r3:492:0.015644444444444443:0.03210550065408781:1,r2:63:0.005777777777777779:0.027595814541259815:1,r1:21:0.002888888888888889:0.03269135509654675:1,r0:14:0.0013333333333333335:0.03347975613277858:1 (r19:3.51759,(((((r18:0.113701,r17:0.113701)0:0.200419,r16:0.31412)0:0.954493,r15:1.26861)0:1.12827,(r14:0.411902,r13:0.411902)0:1.98498)0:0.0985285,(r12:1.52255,(((r11:0.554776,r10:0.554776)0:0.14544,r9:0.700216)0:0.537912,((r8:0.422905,((r7:0.0129533,r6:0.0129533)0:0.373292,r5:0.386245)0:0.0366606)0:0.654669,((r4:0.000254136,r3:0.000254136)0:1.02896,(r2:0.572668,(r1:0.0912342,r0:0.0912342)0:0.481433)0:0.456547)0:0.0483598)0:0.160554)0:0.284423)0:0.972855)0:1.02218);

We provide a convenience function in the API mode for parsing the output file
format (see `iBioGen.util.load_sims()).

## Default CLI args
The default CLI will parse a handful of universally useful arguments:
* `-n`  This is the flag to create a new params file
* `-p`  The flag to specify a params file for a new run of the app
* `-s`  How many simulations to run
* `-c`  How many cores to spin up with the ipyparallel backend
* `-f`  Force the operation (overwrite anything that already exists)
* `-v`  Print out more progress info
* `-q`  Don't print anything to standard out
* `-d`  Turn on debug mode to log debug info to a file
* `-V`  Print version info and exit

Long form arguments:

* `--ipcluster <cluster_id>`    Pass in the cluster ID of a running ipcluster instance

## (Advanced Topics) iBioGen API R bindings
The iBioGen native client (installed through conda) is a standalone command line
program, but it offers a rich API mode which is available by using the R/python
interoperability library [reticluate](https://rstudio.github.io/reticulate/).
The API mode gives much more flexibility and allows for reproducibility in
RMarkdown or juypter.

### Serial Simulations
In the simplest case, you may run simulations serially:

    library(reticulate)

    # import the iBioGen python module
    iBioGen <- import("iBioGen")

    # Create a new iBioGen `Core` object passing in a name
    core = iBioGen$Core("watdo")

    # Print the default parameters
    core$get_params(verbose=TRUE)

    # Set a few parameters
    core$set_param("ClaDS", TRUE)
    core$set_param("ClaDS_alpha", 0.8)

    # In the simplest case you can call the simulate() method, passing in
    # the number of simulations you'd like to perform.
    core$simulate(nsims=2)

    # By default simulations are written to a file, but these can be easily
    # loaded into rstudio with the load_sims() command. load_sims() returns
    # a 'tuple', with the first element being the simulation parameters per
    # sim and the second element being the simulated data.
    res = core$load_sims()

    # Access the simulation parameters
    print(res[1])
