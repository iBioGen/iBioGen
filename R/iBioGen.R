library(reticulate)

# import the iBioGen python module
iBioGen <- import("iBioGen")

simulate <- function(simulation_name = 'my_sim',
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
                     nsims = 1){

    # Create a new iBioGen `Core` object passing in a name
    core = iBioGen$Core("watdo")

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

    res = core$serial_simulate(nsims=nsims)

    return(res)
}
