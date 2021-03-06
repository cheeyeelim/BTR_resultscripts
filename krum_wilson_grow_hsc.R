#(1) Setup paths and environment.
library(BTR)
library(doParallel)

path='~/btr_output/'
setwd(path)

#Setting up for parallel processing.
registerDoParallel(cores=4)

#(2)Load data.
data(krum_bmodel)
data(krum_istate)
data(wilson_raw_data)

bmodel = initialise_model(krum_bmodel)
istate = initialise_data(krum_istate)

#(3)Filtering expression data.
cdata = initialise_raw_data(wilson_raw_data, max_expr = 'low') #do not filter data at this stage. keep the whole matrix.
cell_ind = grepl('cmp', rownames(cdata)) | grepl('gmp', rownames(cdata)) | grepl('mep', rownames(cdata))
fcdata = cdata[cell_ind,] #select only relevant cells.

#(4)Adding extra genes into the model.
extra_genes = setdiff(colnames(wilson_raw_data), bmodel@target)
grown_bmodel = grow_bmodel(extra_genes[c(19, 20)], bmodel) #genes to be added: ldb1, lmo2

#(5)Re-estimate initial state.
#Since CMPs are upsteam of GMPs, and MEPs, the initial state of extra genes can be estimated from CMPs.
tmp_istate = colMeans(cdata[grepl('cmp', rownames(cdata)), extra_genes[c(19, 20)]])
tmp_istate = matrix(round(tmp_istate), nrow=1, dimnames=list(1, names(tmp_istate)))
grown_istate = cbind(istate, tmp_istate)
grown_istate = initialise_data(grown_istate)

result = model_train(cdata=fcdata, bmodel=grown_bmodel, istate=grown_istate, verbose=T)

filename = 'krum_wilson_myeloid_grow_trained_bmodel.rda'
save(result, file=filename)

#Cleaning up.
stopImplicitCluster()
