#(1) Setup paths and environment.
library(BoolTraineR)
library(doParallel)

path='~/bool_final/'
#path='/home/cyl49/ownCloud/Research_Work/working_directory/boolean_project/res1/'
#path='D:/ownCloud/Research_Work/working_directory/boolean_project/res1/'
setwd(path)

inter_bool = T
max_varperrule = 6
acyclic = T
initial_model = T

#Setting up for parallel processing.
registerDoParallel(cores=4) #this automatically calls mclapply() when no cl is given.

#Setting random seed for reproducibility.
set.seed(1)

test_result = c()
for(file_ind in 1:5)
{
  #Load data.
  load(paste(path, ifelse(acyclic, 'acyclic', 'cyclic'), '_data/test_data_yeast', file_ind,'_20n30s10g.rda', sep='')) #test_data.
  
  #Making model inference.
  start_time = proc.time()
  if(initial_model)
  {
    ind_i = sample(seq(1,length(test_data$bm_modmodel)), 1)
    ind_j = 10
    result = model_train(cdata=test_data$cdata, ddata=test_data$ddata, bmodel=test_data$bm_modmodel[[ind_i]][[ind_j]], istate=test_data$istate, and_bool=inter_bool, verbose=T, self_loop=F) #Very time consuming
  } else
  {
    result = model_train(cdata=test_data$cdata, ddata=test_data$ddata, istate=test_data$istate, and_bool=inter_bool, verbose=T, self_loop=F) #Very time consuming
  }
  adj_mat = abs(bm_to_amat(result))
  end_time = proc.time()
  time_taken = end_time - start_time
  
  #Perform model validation.
  true_adj_mat = abs(test_data$true_amat)
  val_res = validate_adjmat(adj_mat, true_adj_mat)
  perf_score = calc_roc(val_res)
  
  #Record results.
  tmp = c(list(true_adj_mat=true_adj_mat), list(adj_mat=adj_mat), list(perf_score=perf_score), list(time_taken=time_taken))
  test_result = c(test_result, list(tmp))
}

#Open file connection for writing output.
file_name = paste('result_btr-', ifelse(initial_model, 'wi', 'wo'), '_train_', ifelse(acyclic, 'acyclic', 'cyclic'), '_yeast.rda', sep='')
save(test_result, file=file_name)

#Cleaning up.
stopImplicitCluster()
