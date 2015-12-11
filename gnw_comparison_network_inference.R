library(BTR)
library(ggplot2)
library(gridExtra) #for plotting multiple graphs
library(colorspace)
library(reshape2)
library(Hmisc) #for capitalisation.

#(1) Load results and combine them.
inter_bool = T
max_varperrule = 6
acyclic = 'cyclic' #acyclic, cyclic, both

# sect_1 ------------------------------------------------------------------

path='~/res1/'
setwd(path)

if(acyclic=='acyclic')
{
  files = list.files(path=path, pattern='^result.+_train_acyclic_yeast[.]rda', full.names = T)
} else if(acyclic=='cyclic')
{
  files = list.files(path=path, pattern='^result.+_train_cyclic_yeast[.]rda', full.names = T)
} else if(acyclic=='both')
{
  files = list.files(path=path, pattern='^result.+_train_[a-z]+_yeast[.]rda', full.names = T)
}

# sect_2 ------------------------------------------------------------------

all_true_adj_mat = list()
all_adj_mat = list()
all_perf_score = c()
all_time_taken = list()
for(i in 1:length(files))
{
  if(acyclic=='acyclic')
  {
    tmp_name = gsub('.+result_(.+)_train_acyclic_yeast[.]rda', '\\1', files[i]) #extract name of method.
  } else if(acyclic=='cyclic')
  {
    tmp_name = gsub('.+result_(.+)_train_cyclic_yeast[.]rda', '\\1', files[i]) #extract name of method.
  } else if(acyclic=='both')
  {
    tmp_name = gsub('.+result_(.+)_train_.+_yeast[.]rda', '\\1', files[i]) #extract name of method.
  }
  
  load(files[i]) #test_result. "true_adj_mat", "adj_mat", "perf_score", "time_taken" 
  
  algo_true_adj_mat = list()
  algo_adj_mat = list()
  algo_perf_score = c()
  algo_time_taken = list()
  for(j in 1:length(test_result))
  {
    algo_true_adj_mat = c(algo_true_adj_mat, list(test_result[[j]]$true_adj_mat))
    algo_adj_mat = c(algo_adj_mat, list(test_result[[j]]$adj_mat))
    algo_perf_score = c(algo_perf_score, list(test_result[[j]]$perf_score))
    algo_time_taken = c(algo_time_taken, list(as.vector(test_result[[j]]$time_taken)))
  }
  algo_perf_score = do.call(rbind, algo_perf_score)
  algo_time_taken = do.call(rbind, algo_time_taken)
  
  if(tmp_name %in% names(all_true_adj_mat)) #only happen in acyclic=both.
  {
    all_true_adj_mat[[tmp_name]] = c(all_true_adj_mat[[tmp_name]], algo_true_adj_mat)
    all_adj_mat[[tmp_name]] = c(all_adj_mat[[tmp_name]], algo_adj_mat)
    all_perf_score[[tmp_name]] = rbind(all_perf_score[[tmp_name]], algo_perf_score)
    all_time_taken[[tmp_name]] = rbind(all_time_taken[[tmp_name]], algo_time_taken)
  } else
  {
    all_true_adj_mat = c(all_true_adj_mat, setNames(list(algo_true_adj_mat), tmp_name))
    all_adj_mat = c(all_adj_mat, setNames(list(algo_adj_mat), tmp_name))
    all_perf_score = c(all_perf_score, setNames(list(algo_perf_score), tmp_name))
    all_time_taken = c(all_time_taken, setNames(list(algo_time_taken), tmp_name))
  }
}

# sect_3 ------------------------------------------------------------------

#Prepare data for plotting.
method_name = toupper(as.vector(sapply(names(all_perf_score), function(x) rep(x, nrow(all_perf_score[[1]])))))
network_name = rep(paste('net_', seq(1,nrow(all_perf_score[[1]])), sep=''), length(all_perf_score))
comb_df = do.call(rbind, all_perf_score)

stopifnot(nrow(comb_df)==length(method_name))

comb_df = data.frame(method_name, network_name, comb_df, stringsAsFactors=F)
comb_df[is.na(comb_df)] = 0
colnames(comb_df) = c(colnames(comb_df)[1:2], toupper(colnames(comb_df)[3:10]))
colnames(comb_df)[10] = 'F-score'

#Setting up consistent colour.
plot_col = rainbow_hcl(length(unique(comb_df[,'method_name'])), alpha=0.8)
names(plot_col) = unique(comb_df[,'method_name'])

plot_df = melt(comb_df)

#1st plot.
comb_p1 = ggplot(plot_df, aes(x=plot_df[,'method_name'], y=plot_df[,"value"])) +
  geom_boxplot() + geom_jitter() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) + 
  xlab('Inference algorithms') + ylab('') + ggtitle(capitalize(acyclic)) +
  facet_wrap(~variable, scales='free_y', ncol=2)

png(paste('all_algorithms_', acyclic, '_score.png', sep=''), width=3000, height=4000, res=300)
plot(comb_p1)
dev.off()

#2nd plot
tmp_list = list()
for(i in unique(method_name))
{
  tmp_row = c(i, colMeans(subset(comb_df[,c(-1,-2)], comb_df[,'method_name']==i)), 
              apply(subset(comb_df[,c(-1,-2)], comb_df[,'method_name']==i), 2, function(x) sd(x)/sqrt(length(x))))
  #tmp_row = c(i, apply(subset(comb_df[,c(-1,-2)], comb_df[,'method_name']==i), 2, sum))
  names(tmp_row) = c('method_name', 'p', 'np', 'r', 'a', 's', 'plr', 'nlr', 'f', 'p_se', 'np_se', 'r_se', 'a_se', 's_se', 'plr_se', 'nlr_se', 'f_se')
  tmp_list = c(tmp_list, list(tmp_row))
}
combsum_df = do.call(rbind, tmp_list)
tmp_rownames = combsum_df[,1]
combsum_df = combsum_df[,-1]
combsum_df = apply(combsum_df, 2, as.numeric)
combsum_df = data.frame(combsum_df)
combsum_df = cbind(method_name=tmp_rownames, combsum_df)

plot_mid_df = melt(combsum_df)
plot_val_df = subset(plot_mid_df, !grepl('[a-z]_se', plot_mid_df[,2]))
plot_se_df = subset(plot_mid_df, grepl('[a-z]_se', plot_mid_df[,2]))
plot_val_df = cbind(plot_val_df, plot_val_df[,3]-plot_se_df[,3], plot_val_df[,3]+plot_se_df[,3])
colnames(plot_val_df) = c('method_name', 'variable', 'value', 'low', 'high')

combmid_p1 = ggplot(plot_val_df, aes(x=plot_val_df[,'method_name'], y=plot_val_df[,'value'], fill=method_name)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=plot_val_df[,"low"], ymax=plot_val_df[,"high"]), width=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) + guides(fill=FALSE) +
  scale_fill_manual("Legend", values = plot_col) +
  xlab('Inference algorithms') + ylab('') + ggtitle(capitalize(acyclic)) +
  facet_wrap(~variable, scales='free_y', ncol=2)

png(paste('all_algorithms_', acyclic, '_meanscore.png', sep=''), width=3000, height=4000, res=300)
plot(combmid_p1)
dev.off()

#Individual plots for F-scores
comb_p8 = ggplot(comb_df, aes(x=reorder(comb_df[,'method_name'], -comb_df[,"F-score"], FUN=median), y=comb_df[,"F-score"])) +
  geom_boxplot() + geom_jitter() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) + 
  xlab('Inference algorithms') + ylab('F-score') + ggtitle(capitalize(acyclic))

png(paste('all_algorithms_', acyclic, '_fscore.png', sep=''), width=2000, height=2000, res=300)
plot(comb_p8)
dev.off()

combmid_p8 = ggplot(combsum_df, aes(x=reorder(combsum_df[,'method_name'], -combsum_df[,"f"]), y=f, fill=method_name)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=combsum_df[,"f"] - combsum_df[,'f_se'], ymax=combsum_df[,"f"] + combsum_df[,'f_se']), width=.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) + guides(fill=FALSE) +
  scale_fill_manual("Legend", values = plot_col) +
  xlab('Inference algorithms') + ylab('F-score') + ggtitle(capitalize(acyclic))

png(paste('all_algorithms_', acyclic, '_fmeanscore.png', sep=''), width=2000, height=2000, res=300)
plot(combmid_p8)
dev.off()
