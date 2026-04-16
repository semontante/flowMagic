# Remember to load the libraries to read, visualize and analyse flow cytometry data, including flowMagic.
library(flowMagic)
library(flowCore)
library(ggplot2)
library(flowWorkspace)
library(devtools)
library(CytoML)

# Import data and preprocessing ####

# Import FCS files from flowMagic extdata directory
path_dir <- system.file("extdata/fcs_files", package = "flowMagic")
fs <- read.flowSet(
  path = path_dir,
  transformation = FALSE,
  pattern = ".fcs"
)
# Create a GatingSet
gs <- GatingSet(fs)
# View sample names
sampleNames(gs)
# Load compensation matrix from extdata
path_comp_file <- system.file(
  "extdata/comp_matrix.csv",
  package = "flowMagic"
)
comp_matrix <- read.csv(path_comp_file, check.names = FALSE)
# Apply compensation to the GatingSet
gs <- compensate(gs, comp_matrix)
# Channels to be transformed using a logicle transform
channels_to_transform <- c(
  "FL1.A", "FL2.A", "FL3.A", "Viability_dye"
)
# Estimate transform based on the first sample
trans_list <- estimateLogicle(
  gs[[1]],
  channels = channels_to_transform
)
# Apply transformation
gs <- transform(gs, trans_list)
# visualize pre-processed data
flowMagic::magicPlot_fs(
  fs = fs,
  sample_id = 1,
  channel_x = "FSC.A",
  channel_y = "FSC.H"
)
# manual gate singlets ####
list_out<-magicGating(fs = fs,sample_id = c("sim_standard_gating_01.fcs"),
                      channel_x = "FSC.A",channel_y = "FSC.H",size_points=0.5)

df_1<-list_out$list_gated_data$sim_standard_gating_01.fcs
magicPlot(df = df_1) # visually check manually gated plot
# Apply the same polygonGate to all samples in the GatingSet
gs_pop_add(gs, list_out$list_poly_gates$sim_standard_gating_01.fcs, parent = "root", name = "Singlets")
# Recompute the gating
recompute(gs)
# plot gating tree
plot(gs)
# we can save the new GatingSet with the new gating information
# Note that the path to an empty directory is required.
# In this example, the exported_gs_directory is a new empty directory.
save_gs(gs = gs,path = "/home/rstudio/main/Data/sim_data/exported_gs")


# Automated gating live cell using Template model ####
# import gs (can also be imported from flowJo or FCS express)
gs<-load_gs(path = "/home/rstudio/main/Data/sim_data/exported_gs")
sampleNames(gs)

ff<-flowMagic::get_flowframe_from_gs(gs = gs,node_name = "Singlets",sample_id = "sim_standard_gating_02.fcs")
flowMagic::magicPlot_fs(fs = ff,channel_x = "Viability_dye",channel_y = "FSC.A")

ff<-flowMagic::get_flowframe_from_gs(gs = gs,node_name = "Singlets",sample_id = "sim_standard_gating_05.fcs")
flowMagic::magicPlot_fs(fs = ff,channel_x = "Viability_dye",channel_y = "FSC.A")

# make templates
list_out_1<-magicGating(fs = gs,sample_id = c("sim_standard_gating_02.fcs","sim_standard_gating_05.fcs"),
                        channel_x = "Viability_dye",channel_y = "FSC.A",gs_node = "Singlets",label_pol = "1")

list_out_2<-magicGating(fs = gs,sample_id = c("sim_standard_gating_02.fcs","sim_standard_gating_05.fcs"),
                        channel_x = "Viability_dye",channel_y = "FSC.A",gs_node = "Singlets",label_pol = "2")

list_out_final<-merge_magicGating_labels(list_out_1 = list_out_1,list_out_2 = list_out_2)

# we convert these gated samples to a training set
ref_train<-get_train_data(paths_file = list_out_final)
# training step based on template data to generate template model
ref_model_info<-magicTrain(df_train = ref_train,train_model = "rf")

# get test data
list_test_data<-flowMagic::export_raw_gs_plots(gs = gs,node_name = "Singlets",channel_x = "Viability_dye",channel_y
                                               = "FSC.A",
                                               return_data = T)

# perform automated gating based on template model
list_dfs_pred<-magicPred_all(list_test_data = list_test_data,magic_model = NULL,ref_data_train = ref_train,
                             ref_model_info = ref_model_info,n_cores = 1)
# export gated plots
exports_plots(list_gated_data = list_dfs_pred,path_output = "~/main/Data/results_sim_data")
# visualize all gated plots wrapped together
flowMagic::magic_plot_wrap(list_gated_data = list_dfs_pred,n_col_wrap = 3,
                           size_points=0.5,size_title_x=15,size_title_y=15,size_axis_text=10)

# try changing the gates based on the type of data in the template (e.g., changing gates of high density vs low density populations)

# add gated live cells into Gating hierarchy
list_poly_gates<-flowMagic::flowmagic_pred_to_poly_gates(list_df = list_dfs_pred,gate_label = "Live_cells",pred_label = "1")
gs<-flowMagic::flowmagic_pred_to_gs(list_poly_gates = list_poly_gates,gs = gs,parent_node = "Singlets")
list_poly_gates<-flowmagic_pred_to_poly_gates(list_df = list_dfs_pred,gate_label = "Dead_cells",pred_label = "2")
gs<-flowMagic::flowmagic_pred_to_gs(list_poly_gates = list_poly_gates,gs = gs,parent_node = "Singlets")
save_gs(gs = gs,path = "/home/rstudio/main/Data/sim_data/exported_gs_2")

# visualize hierarchy with flowMagic
out<-get_hierarchy_all_pops(gh=gs[[1]],export_visnet = F)

out$visnet

# plot gating tree using flowCore/flowWorkspace
plot(gs)

# in flowMagic hierarchical tree, populations defined by same pair of markers have the same color.

# Automated gating live cells using Generalized model ####
# we need to reload the previous GatingSet we did not modify with the new Live cells gate
gs<-load_gs(path = "/home/rstudio/main/Data/sim_data/exported_gs")

# First, we need to load the generalized model (divided in model A and B) from the Federated Research Data Repository (FRDR).
# Once downloaded we need to load them in R using the readRDS() function.
model_a<-readRDS("~/main/GP_model/models_trained_n_gates_final/models_trained_to_predict_n_gates_final/training_rf_index_3000train10val_2ntree_500points_100folds_31000_consensus_plots_pred_n_gates.RData")

model_b<-readRDS("~/main/GP_model/models_trained_n_gates_final/models_trained_to_predict_classes_final/list_models_all_n_gates.RData")

# get test data
list_test_data<-flowMagic::export_raw_gs_plots(gs = gs,node_name = "Singlets",channel_x = "Viability_dye",channel_y
                                               = "FSC.A",
                                               return_data = T)

# perform automated gating based on generalized model

# with a single predefined number of gates for all samples
list_dfs_pred<-magicPred_all(list_test_data = list_test_data,magic_model = model_b,magic_model_n_gates = 2,
                             n_cores = 8,n_points_per_plot = 15000)

# with a single predefined number of gates for selected samples
list_dfs_pred<-magicPred_all(list_test_data = list_test_data,sample_id = c(1,5,10),magic_model = model_b,magic_model_n_gates = 2,
                             n_cores = 8,n_points_per_plot = 15000) # instead of positional indices, samples names can also be indicated

# with different predefined number of gates based on the samples
# generating dataframe that pairs sample names with predefined number of gates.
n_gates_samples<-rep(2,length(list_test_data))
sample_name<-names(list_test_data)
n_gates_df<-as.data.frame(cbind(sample_name,n_gates_samples))
n_gates_df$n_gates_samples[5]<-1
n_gates_df$n_gates_samples[6]<-1
n_gates_df$n_gates_samples[9]<-1
list_dfs_pred<-magicPred_all(list_test_data = list_test_data,magic_model = model_b,
                             n_gates_df=n_gates_df,
                             n_cores = 8,n_points_per_plot = 15000)

# with number of gates predicted by model A for samples
list_dfs_pred<-magicPred_all(list_test_data = list_test_data,magic_model = model_b,magic_model_n_gates = model_a,
                             n_cores = 8,n_points_per_plot = 15000,thr_dist = 0.2)

# visualize all gated plots wrapped together
flowMagic::magic_plot_wrap(list_gated_data = list_dfs_pred,n_col_wrap = 3,
                           size_points=0.5,size_title_x=15,size_title_y=15,size_axis_text=10)

magicPlot(list_dfs_pred$sim_standard_gating_02.fcs$df_test_original)

# visualize density of events next to each axis.
magicPlot(list_dfs_pred$sim_standard_gating_05.fcs$df_test_original,show_marginals = T)
magicPlot(list_dfs_pred$sim_standard_gating_02.fcs$df_test_original,show_marginals = T)
# export gated plots
exports_plots(list_gated_data = list_dfs_pred,path_output = "~/main/Data/results_sim_data")
