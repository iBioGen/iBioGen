
#' Simulate a phylogeny with abundance and genetic diversity at the tips
#'
#' @details
#' A birth/death process with abundances and genetic diversity at the tips.
#' Abundance can evolve either as a BM process with random fission at speciation
#' events, or the rate of change (r) of abundance can evolve as BM in which case
#' abundance (n) changes through time (dt) via n * exp(r*dt). Speciation rate
#' can also shift at branching events in the manner of ClaDS.
#'
#' The iBioGen package is a product of the iBioGen project, an EU-funded
#' (Horizon 2020) project, focused on island biodiversity genomics. More details
#' about the project can be found on the website: https://www.ibiogen.eu/
#' @param simulation_name The name of this simulation scenario
#' @param project_dir Where to save files.
#' @param birth_rate Speciation rate
#' @param stop_criterion Whether to stop on ntaxa or time
#' @param ntaxa Number of taxa to simulate if sto is `ntaxa`
#' @param time Amount of time to simulate if stop is `time`
#' @param process Whether to evolve `abundance` or growth `rate` via BM
#' @param ClaDS Whether to allow speciation rates to change along the branches a la ClaDS
#' @param abundance_mean Ancestral abundance at time 0
#' @param abundance_sigma Rate at which abundance changes if process is `abundance`
#' @param growth_rate_mean Ancestral population growth rate at time 0
#' @param growth_rate_sigma Rate at which growth rate changes if process is `rate`
#' @param ClaDS_sigma Rate at which speciation rate changes if ClaDS is True
#' @param ClaDS_alpha Rate shift if ClaDS is True
#' @param sequence_length Length of the genomic region simulated, in base pairs
#' @param mutation_rate Mutation rate per base per generation
#' @param sample_size Number of samples to draw for calculating genetic diversity
#' @param abundance_scaling Scaling abundance to Ne. Can be None, log, ln or a ratio
#' @param nsims Number of independent simulations to perform
#' @param quiet Suppress printing the progress bar
#' @export
#' @examples
#' sim_tree() 
#' sim_tree(ClaDS=TRUE, ntaxa=200)
sim_tree <- function(simulation_name = 'my_sim',
                     project_dir = './default_iBioGen',
                     birth_rate = 1,
                     stop_criterion = "taxa",
                     ntaxa = 20,
                     time = 4,
                     process = "rate",
                     ClaDS = FALSE,
                     abundance_mean = 500000,
                     abundance_sigma = 0.1,
                     growth_rate_mean = 0,
                     growth_rate_sigma = 0.01,
                     ClaDS_sigma = 0.9,
                     ClaDS_alpha = 0.1,
                     sequence_length = 1000,
                     mutation_rate = 0.0000001,
                     sample_size = 10,
                     abundance_scaling = "None",
                     nsims = 1,
                     quiet = FALSE){

    library(reticulate)
    tryCatch({
        conda_version()
    }, error = function(err) {
        print("  First run configuring iBioGen library.")
        print("  This could take 5-10 minutes.")
        print("  Installing miniconda")
        install_miniconda()
        print("  Installing iBioGen backend.")
        conda_install("r-reticulate",
                      packages="ibiogen",
                      channel=c("conda-forge", "iBioGen"),
                      python_version=3.7)
    })

    # import the iBioGen python module
    tryCatch({
        pyBioGen <- import("iBioGen")
    }, error = function(err) {
        print("  Error in loading iBioGen backend.")
        print(paste(" ", err))
        # Do sttuff here
    })

    # Create a new iBioGen `Core` object passing in a name
    core = pyBioGen$Core("watdo")

    # Set parameters
    core$set_param('simulation_name', simulation_name)
    core$set_param('project_dir', project_dir)
    core$set_param('birth_rate', birth_rate)
    core$set_param('stop_criterion', stop_criterion)
    core$set_param('ntaxa', ntaxa)
    core$set_param('time', time)
    core$set_param('process', process)
    core$set_param('ClaDS', ClaDS)
    core$set_param('abundance_mean', abundance_mean)
    core$set_param('abundance_sigma', abundance_sigma)
    core$set_param('growth_rate_mean', growth_rate_mean)
    core$set_param('growth_rate_sigma', growth_rate_sigma)
    core$set_param('ClaDS_alpha', ClaDS_alpha)
    core$set_param('ClaDS_sigma', ClaDS_sigma)
    core$set_param('sequence_length', sequence_length)
    core$set_param('mutation_rate', mutation_rate)
    core$set_param('sample_size', sample_size)
    core$set_param('abundance_scaling', abundance_scaling)

    res = core$serial_simulate(nsims=nsims, quiet=quiet)

    return(res)
}
