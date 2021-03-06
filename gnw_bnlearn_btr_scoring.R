#(1) Setup paths and environment.
#library(BoolTraineR)
library(BTR)
library(bnlearn)
library(ggplot2)
library(gridExtra) #for plotting multiple graphs
library(reshape2)

path='~/btr_output/'
setwd(path)

max_varperrule = 6
nonoise=F
subset_bool = F

for(acyclic in c(T,F))
{
  for(single_cell in c(T,F))
  {
    for(file_ind in 1:5)
    {
      #Load data.
      load(paste(path, ifelse(acyclic, 'acyclic', 'cyclic'), '_data/test_data_', ifelse(nonoise, 'nonoise_', ''),'yeast', file_ind,'_20n30s10g.rda', sep=''))  #test_data.
      
      #score the true model under Bayesian framework.
      if(acyclic)
      {
        true_ndamat = abs(test_data$true_amat)
        true_bngraph = empty.graph(colnames(test_data$true_amat))
        amat(true_bngraph) = true_ndamat #if this step throws error, that means the true graph is not DAG. Has to regenerate the random data using GNW.
        true_bnscore = score(true_bngraph, as.data.frame(test_data$cdata))
      }
      
      #score the true model under Boolean model framework.
      true_bmmodel = amat_to_bm(test_data$true_amat)
      overlap_gene = unname(colnames(test_data$cdata))
      if(single_cell)
      {
        if(subset_bool)
        {
          set.seed(1)
          sind = sample(1:nrow(test_data$sc_cdata), 50)
          fcdata = test_data$sc_cdata[sind, ]
        } else
        {
          fcdata = test_data$sc_cdata
        }
      } else
      {
        if(subset_bool)
        {
          set.seed(1)
          sind = sample(1:nrow(test_data$cdata), 50)
          fcdata = test_data$cdata[sind, ]
        } else
        {
          fcdata = test_data$cdata
        }
      }
      true_bmscore = calc_mscore(bmodel=true_bmmodel, istate=test_data$istate, fcdata=fcdata, overlap_gene=overlap_gene, max_varperrule=max_varperrule)
      
      #score the other bn models under Bayesian framework using GNW data.
      bnmod_gnw_bnscore = c()
      if(acyclic)
      {
        for(i in 1:length(test_data$bn_modamat))
        {
          for(j in 1:length(test_data$bn_modamat[[i]]))
          {
            tmp_bngraph = empty.graph(colnames(test_data$bn_modamat[[i]][[j]]))
            amat(tmp_bngraph) = test_data$bn_modamat[[i]][[j]]
            tmp_score = score(tmp_bngraph, as.data.frame(test_data$cdata))
            #tmp_score = score(tmp_bngraph, as.data.frame(bm_cdata))
            names(tmp_score) = test_data$bn_step[[i]][[j]]
            
            bnmod_gnw_bnscore = c(bnmod_gnw_bnscore, tmp_score)
          }
        }
        bnmod_gnw_bnscore = split(bnmod_gnw_bnscore, names(bnmod_gnw_bnscore))
        bnmod_gnw_bnscore = bnmod_gnw_bnscore[order(as.numeric(names(bnmod_gnw_bnscore)))]
        bnmod_gnw_bnscore = do.call(cbind, bnmod_gnw_bnscore)
        #Add in true score.
        bnmod_gnw_bnscore = cbind('0'=rep(true_bnscore, nrow(bnmod_gnw_bnscore)), bnmod_gnw_bnscore)
        rownames(bnmod_gnw_bnscore) = NULL
        colnames(bnmod_gnw_bnscore) = paste(seq(0,40,2))
      }
      
      #score the other bm models under Boolean model framework using GNW data.
      bmmod_gnw_bmscore = c()
      y_bmscore = c()
      za_bmscore = c()
      zb_bmscore = c()
      for(i in 1:length(test_data$bm_modmodel))
      {
        for(j in 1:length(test_data$bm_modmodel[[i]]))
        {
          
          tmp_bmmodel = test_data$bm_modmodel[[i]][[j]]
          overlap_gene = unname(colnames(test_data$cdata))
          tmp_bmmodel@target = overlap_gene
          tmp_score = calc_mscore(bmodel=tmp_bmmodel, istate=test_data$istate, fcdata=fcdata, overlap_gene=overlap_gene, max_varperrule=max_varperrule, detail=T)
          
          names(tmp_score) = rep(j, length(tmp_score))
          
          #c('f', 'y', 'zb')
          bmmod_gnw_bmscore = c(bmmod_gnw_bmscore, tmp_score[1])
          y_bmscore = c(y_bmscore, tmp_score[2])
          za_bmscore = c(za_bmscore, tmp_score[3])
          zb_bmscore = c(zb_bmscore, tmp_score[4])
        }
      }
      bmmod_gnw_bmscore = split(bmmod_gnw_bmscore, names(bmmod_gnw_bmscore))
      bmmod_gnw_bmscore = bmmod_gnw_bmscore[order(as.numeric(names(bmmod_gnw_bmscore)))]
      bmmod_gnw_bmscore = do.call(cbind, bmmod_gnw_bmscore)
      #Add in true score.
      bmmod_gnw_bmscore = cbind('0'=rep(true_bmscore, nrow(bmmod_gnw_bmscore)), bmmod_gnw_bmscore)
      rownames(bmmod_gnw_bmscore) = NULL
      colnames(bmmod_gnw_bmscore) = paste(seq(0,40,2))
      
      y_bmscore = split(y_bmscore, names(y_bmscore))
      y_bmscore = y_bmscore[order(as.numeric(names(y_bmscore)))]
      y_bmscore = do.call(cbind, y_bmscore)
      rownames(y_bmscore) = NULL
      colnames(y_bmscore) = paste(seq(2,40,2))
      
      za_bmscore = split(za_bmscore, names(za_bmscore))
      za_bmscore = za_bmscore[order(as.numeric(names(za_bmscore)))]
      za_bmscore = do.call(cbind, za_bmscore)
      rownames(za_bmscore) = NULL
      colnames(za_bmscore) = paste(seq(2,40,2))
      
      zb_bmscore = split(zb_bmscore, names(zb_bmscore))
      zb_bmscore = zb_bmscore[order(as.numeric(names(zb_bmscore)))]
      zb_bmscore = do.call(cbind, zb_bmscore)
      rownames(zb_bmscore) = NULL
      colnames(zb_bmscore) = paste(seq(2,40,2))
      
      save(bnmod_gnw_bnscore, bmmod_gnw_bmscore, y_bmscore, za_bmscore, zb_bmscore, file=paste(ifelse(acyclic, 'acyclic', 'cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_yeast', file_ind, '_scores.rda', sep=''))
      #save(bnmod_gnw_bnscore, bmmod_gnw_bmscore, y_bmscore, zb_bmscore, file=paste(ifelse(acyclic, 'acyclic', 'cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_yeast', file_ind, '_scores.rda', sep=''))
    }
    
    ##################################################################################################################################
    
    #Setup the true networks.
    true_model_list = list()
    for(file_ind in 1:5)
    {
      #Load data.
      load(paste(path, ifelse(acyclic, 'acyclic', 'cyclic'), '_data/test_data_', ifelse(nonoise, 'nonoise_', ''),'yeast', file_ind,'_20n30s10g.rda', sep=''))  #test_data.
      
      true_bmmodel = amat_to_bm(test_data$true_amat)
      true_model_list = c(true_model_list, list(true_bmmodel))
    }
    
    #Make plots of true models.
    for(i in 1:length(true_model_list))
    {
      png(paste(i, ifelse(acyclic, '_acyclic', '_cyclic'), '_true_network.png', sep=''), width=2000, height=2000, res=300)
      plotBM(true_model_list[[i]])
      title(paste('True ', ifelse(acyclic, 'acyclic', 'cyclic'), ' network', i), cex.main=3)
      dev.off()
    }
    
    #Setup the scores for modified models.
    bnscore_list = list()
    bmscore_list = list()
    for(file_ind in 1:5)
    {
      load(paste(ifelse(acyclic, 'acyclic', 'cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_yeast', file_ind, '_scores.rda', sep='')) #bnmod_gnw_bnscore, bmmod_gnw_bmscore

      bnscore_list = c(bnscore_list, list(bnmod_gnw_bnscore))
      bmscore_list = c(bmscore_list, list(bmmod_gnw_bmscore))
    }
    
    if(acyclic)
    {
      bnscore_df = melt(bnscore_list)
      bnscore_df = bnscore_df[,-1]
      bnscore_df[,3] = paste('Network', bnscore_df[,3])
      colnames(bnscore_df) = c('steps', 'score', 'network')
      
      bnscore_mid_df = melt(lapply(bnscore_list, function(x) colMeans(x)))
      bnscore_mid_df[,2] = paste('Network', bnscore_mid_df[,2])
      bnscore_mid_df = cbind(rep(seq(0,40,2), 5), bnscore_mid_df)
      
      bnscore_se = melt(lapply(bnscore_list, function(x) apply(x, 2, function(x) sd(x)/sqrt(length(x)))))
      bnscore_se = bnscore_se[,1,drop=F]
      bnscore_mid_df = cbind(bnscore_mid_df, bnscore_mid_df[,2]-bnscore_se[,1], bnscore_mid_df[,2]+bnscore_se[,1])
      colnames(bnscore_mid_df) = c('steps', 'score', 'network', 'low', 'high')
    }

    bmscore_df = melt(bmscore_list)
    bmscore_df = bmscore_df[,-1]
    bmscore_df[,3] = paste('Network', bmscore_df[,3])
    colnames(bmscore_df) = c('steps', 'score', 'network')
    
    bmscore_mid_df = melt(lapply(bmscore_list, function(x) colMeans(x)))
    bmscore_mid_df[,2] = paste('Network', bmscore_mid_df[,2])
    bmscore_mid_df = cbind(rep(seq(0,40,2), 5), bmscore_mid_df)
    
    bmscore_se = melt(lapply(bmscore_list, function(x) apply(x, 2, function(x) sd(x)/sqrt(length(x)))))
    bmscore_se = bmscore_se[,1,drop=F]
    bmscore_mid_df = cbind(bmscore_mid_df, bmscore_mid_df[,2]-bmscore_se[,1], bmscore_mid_df[,2]+bmscore_se[,1])
    colnames(bmscore_mid_df) = c('steps', 'score', 'network', 'low', 'high')
    
    #Make plot objects of scoring functions.
    if(acyclic)
    {
      p1_bn_box = ggplot(bnscore_df, aes(x=factor(bnscore_df[,'steps']), y=bnscore_df[,'score'])) +
        geom_boxplot() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BIC scoring function') +
        scale_x_discrete(labels=unique(bnscore_df[,'steps'])) +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20), axis.text.x = element_text(size=10))
      
      p2_bm_box = ggplot(bmscore_df, aes(x=factor(bmscore_df[,'steps']), y=bmscore_df[,'score'])) +
        geom_boxplot() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BSS scoring function') +
        scale_x_discrete(labels=unique(bmscore_df[,'steps'])) +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20), axis.text.x = element_text(size=10))
      
      png(paste('boolbaye', ifelse(acyclic, '_acyclic', '_cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_boxplot_compare_score.png', sep=''), width=3000, height=5000, res=300)
      grid.arrange(p1_bn_box, p2_bm_box, ncol=2)
      dev.off()
      
      p1_bn_mid = ggplot(bnscore_mid_df, aes(x=bnscore_mid_df[, 'steps'], y=bnscore_mid_df[,'score'])) +
        geom_errorbar(aes(ymin=bnscore_mid_df[, 'low'], ymax=bnscore_mid_df[, 'high'])) +
        geom_line() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BIC scoring function') +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20))
      
      p2_bm_mid = ggplot(bmscore_mid_df, aes(x=bmscore_mid_df[, 'steps'], y=bmscore_mid_df[,'score'])) +
        geom_errorbar(aes(ymin=bmscore_mid_df[, 'low'], ymax=bmscore_mid_df[, 'high'])) +
        geom_line() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BSS scoring function') +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20))
      
      png(paste('boolbaye', ifelse(acyclic, '_acyclic', '_cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_mean_compare_score.png', sep=''), width=3000, height=5000, res=300)
      grid.arrange(p1_bn_mid, p2_bm_mid, ncol=2)
      dev.off()
    } else
    {
      p2_bm_box = ggplot(bmscore_df, aes(x=factor(bmscore_df[,'steps']), y=bmscore_df[,'score'])) +
        geom_boxplot() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BSS scoring function') +
        scale_x_discrete(labels=unique(bmscore_df[,'steps'])) +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20), axis.text.x = element_text(size=10))
      
      png(paste('boolbaye', ifelse(acyclic, '_acyclic', '_cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_boxplot_compare_score.png', sep=''), width=3000, height=5000, res=300)
      grid.arrange(p2_bm_box, ncol=2)
      dev.off()
      
      p2_bm_mid = ggplot(bmscore_mid_df, aes(x=bmscore_mid_df[, 'steps'], y=bmscore_mid_df[,'score'])) +
        geom_errorbar(aes(ymin=bmscore_mid_df[, 'low'], ymax=bmscore_mid_df[, 'high'])) +
        geom_line() + xlab('Number of different edges') + ylab('Scores') + ggtitle('BSS scoring function') +
        facet_wrap(~network, scales='free_y', ncol=1) + 
        theme(text = element_text(size=20))
      
      png(paste('boolbaye', ifelse(acyclic, '_acyclic', '_cyclic'), ifelse(single_cell, '_sc', '_nonsc'), '_mean_compare_score.png', sep=''), width=3000, height=5000, res=300)
      grid.arrange(p2_bm_mid, ncol=2)
      dev.off()
    }
  }
}
