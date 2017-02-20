### phylogeographic model testing with mtDNA and nuDNA S. torotoro data

#setwd("/Users/ethanlinck/Dropbox/Syma/syma_phylogeography/")
#install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com");install.packages("snow");install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com");install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")
library(snow); library(BioGeoBEARS); library(optimx);library(gdata);library(parallel)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully

# prune tree for analysis
tree <- read.nexus("timetree.tre")
tree <- drop.tip(tree, "h_sanctus", trim.internal = TRUE)
old_tipnames <- tree$tip.label
new_tipnames <- c("NNG", "ORO", "NBI", "NNG", "NNG", "SNG", "SNG", "SNG", "SNG", "NNG", "WNG", "WPI", "WNG", "WNG", "NNG", "WNG", "WNG", "NNG")
tmpmat <- cbind(old_tipnames, new_tipnames)
df <- as.data.frame(tmpmat, stringsAsfactors=FALSE)
names(df) <- c("OTUs", "species")
pdffn <- "pruning_PDF.pdf"
pdf(pdffn, height=6, width=6)
result <- prune_specimens_to_species(original_tr=tree, xls=df, group_name="default", titletxt="", areas_abbr=NULL, plot_intermediate=TRUE)
dev.off()
cmstr <- paste0("open ",pdffn)
system(cmstr)
pruned_tree <- result$tr
plot(pruned_tree)
write.tree(pruned_tree, file="pruned_tree.newick")
trfn = "pruned_tree.newick"
plot(pruned_tree)
axisPhylo()
title("Divergence of major Syma torotoro mtDNA lineages")

#input geography file
geogfn = "syma_geog.data"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size = 3
numstates_from_numareas(numareas=2, maxareas=2, include_null_range=TRUE)

#input distance file
distfn = "distances.txt"

#Run a DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$distance = distfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001 
BioGeoBEARS_run_object$include_null_range = TRUE 
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = TRUE     
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE 
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "biogeobears.Rdata"
if (runslow)
  {
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
  } else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#PLOT DEC
analysis_titletxt <- "BioGeoBEARS DEC Syma torotoro"
results_object <- res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=pruned_tree, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=pruned_tree, tipranges=tipranges)

