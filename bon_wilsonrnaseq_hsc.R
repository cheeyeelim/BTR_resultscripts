#(1) Setup paths and environment.
library(BTR)
library(doParallel)

path='~/bool_final/'
setwd(path)

#Setting up for parallel processing.
registerDoParallel(cores=4)

#(2)Load data.
data(bon_bmodel)
data(bon_istate)
data(wilson_raw_rnaseq)

inter_bool = T
max_varperrule = 6

bmodel = initialise_model(bon_bmodel)
istate = initialise_data(bon_istate)

#(3)Filtering expression data.
tmp_data = initialise_raw_data(wilson_raw_rnaseq) #do not filter data at this stage. keep the whole matrix.
cdata = tmp_data[[1]] #continuous data
ddata = tmp_data[[2]] #discretised data
fcdata = cdata[grepl('hsc', rownames(cdata)),] #select only HSC cells.
fddata = ddata[grepl('hsc', rownames(ddata)),]

#(3)Filtering expression data.
result = model_train(cdata=fcdata, ddata=fddata, bmodel=bmodel, istate=istate, verbose=T)

filename = 'bon_wilson_rnaseq_HSC_trained_bmodel.rda'
save(result, file=filename)

#Cleaning up.
stopImplicitCluster()
