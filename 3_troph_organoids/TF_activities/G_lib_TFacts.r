library(dorothea)
library(reshape2)


# Generate dorothea regulons object
viper_gset = get(load('utils/dorotheav2-top10scoring_VentoLab20201111.rdata'))
for (i in names(viper_gset))
    viper_gset[[i]] = data.frame(tf=i, confidence=strsplit(i,split = '_')[[1]][2], target=names(viper_gset[[i]]$tfmode), mor=viper_gset[[i]]$tfmode)
dorothea_regulon = melt(viper_gset, id.vars = names(viper_gset[[1]]))[, -5]


# # Generate pathways gsets object
# viper_gset = get(load('~/farm/gsea/genesets/GObpUniprot_viper.rdata'))
# viper_gset = append(x = viper_gset,  get(load('~/farm/gsea/genesets/UniprotKeywords_viper.rdata')))
# viper_gset = append(x = viper_gset,  get(load('~/farm/gsea/genesets/ovarianDevelopment_viper.rdata')))
# for (i in names(viper_gset))
#     viper_gset[[i]] = data.frame(tf=i, confidence='high', target=names(viper_gset[[i]]$tfmode), mor=1)
# pathways_gsa = melt(viper_gset, id.vars = names(viper_gset[[1]]))[, -5]

plot_TFactivities = function(df_TFact, anndataO_doro, outfile = 'heatmap_TFactivities.pdf', TFs=NULL){
            #  Plot TF activities and expression
        require(dplyr)
        require(tibble)
        require(pheatmap)
        require(tidyr)


        # Plot significant TFs
    if (is.null(TFs) ) 
        TFs = subset(df_TFact, avg_logFC > 0 & p_val_adj < 0.05) %>%
          group_by(cluster) %>%
          group_map(~ head(.x, 5L)$gene) %>%
          unlist(.)
        

        # We transform Viper scores, scaled by seurat, into a data frame to better 
        ## handling the results
        viper_scores_df <- GetAssayData(anndataO_doro, slot = "scale.data", 
                                            assay = "dorothea") %>%
          data.frame(check.names = F) %>%
          t()

        viper_scores_df = viper_scores_df[ , unique(TFs)]

        ## We create a data frame containing the cells and their clusters
        CellsClusters <- data.frame(cell = names(Idents(anndataO_doro)), 
                                    cell_type = as.character(Idents(anndataO_doro)),
                                    check.names = F)

        ## We create a data frame with the Viper score per cell and its clusters
        viper_scores_clusters <- viper_scores_df  %>%
          data.frame() %>% 
          rownames_to_column("cell") %>%
          gather(tf, activity, -cell) %>%
          inner_join(CellsClusters)

        ## We summarize the Viper scores by cellpopulation
        summarized_viper_scores <- viper_scores_clusters %>% 
          group_by(tf, cell_type) %>%
          summarise(avg = mean(activity),
                    std = sd(activity))


        ## We prepare the data for the plot
        summarized_viper_scores_df <- summarized_viper_scores %>%
          dplyr::select(-std) %>%   
          spread(tf, avg) %>%
          data.frame(row.names = 1, check.names = FALSE) 

        summarized_viper_scores_df = summarized_viper_scores_df[, gsub('-', '.', unique(TFs))]
        summarized_viper_scores_df = summarized_viper_scores_df[ levels(Idents(anndataO_doro)), ]

        palette_length = 100
        my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

        my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                           length.out=ceiling(palette_length/2) + 1),
                       seq(max(summarized_viper_scores_df)/palette_length, 
                           max(summarized_viper_scores_df), 
                           length.out=floor(palette_length/2)))

        viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                               fontsize_row = 10, 
                               color=my_color, breaks = my_breaks, 
                               main = "TFs", angle_col = 45,
                               treeheight_col = 0,  border_color = NA) 


        # # Add is TF a DE
        # dfDEGs$value = ''
        # dfDEGs$value[ dfDEGs$p_val_adj < 0.05 & dfDEGs$avg_logFC > 0] = '+'
        # dfDEGs$value[ dfDEGs$p_val_adj < 0.05 & dfDEGs$avg_logFC < 0] = '-'
        # m = acast(subset(dfDEGs, gene %in% df_TFact$TF), cluster~gene, fill = '')

        labels = summarized_viper_scores_df
        labels[] = ''
        # tfs = gsub('.AA$', '', colnames(summarized_viper_scores_df)) %>%  gsub('.[A-E]$', '', .)
        # for (tf in tfs){
        #     if (tf %in% colnames(m))
        #         labels[ rownames(m), tfs == tf ] = m[, tf]
        # }
        message('saving csv as:', gsub('.pdf', '.csv', outfile))
        write.csv(summarized_viper_scores_df, file = gsub('.pdf', '.csv', outfile))
        message('saving pdf as:', outfile)
        viper_hmap <- pheatmap(summarized_viper_scores_df, fontsize=14, display_numbers = labels, cluster_rows = F, scale = 'none',
                               cluster_cols = F,
                               cellwidth = 12, cellheight = 10, 
                               filename = outfile,
                               fontsize_row = 10, 
                               color=my_color, breaks = my_breaks, 
                               main = "TFs", angle_col = 90,
                               treeheight_col = 0,  border_color = NA) 

}


get_TFact_cluster = function(anndataO_doro, outfile = 'TFactivities.csv'){
            #  Plot TF activities and expression
        require(dplyr)
        require(tibble)
        require(pheatmap)
        require(tidyr)


        

        # We transform Viper scores, scaled by seurat, into a data frame to better 
        ## handling the results
        viper_scores_df <- GetAssayData(anndataO_doro, slot = "scale.data", 
                                            assay = "dorothea") %>%
          data.frame(check.names = F) %>%
          t()


        ## We create a data frame containing the cells and their clusters
        CellsClusters <- data.frame(cell = names(Idents(anndataO_doro)), 
                                    cell_type = as.character(Idents(anndataO_doro)),
                                    check.names = F)

        ## We create a data frame with the Viper score per cell and its clusters
        viper_scores_clusters <- viper_scores_df  %>%
          data.frame() %>% 
          rownames_to_column("cell") %>%
          gather(tf, activity, -cell) %>%
          inner_join(CellsClusters)

        ## We summarize the Viper scores by cellpopulation
        summarized_viper_scores <- viper_scores_clusters %>% 
          group_by(tf, cell_type) %>%
          summarise(avg = mean(activity),
                    std = sd(activity))


        ## We prepare the data for the plot
        summarized_viper_scores_df <- summarized_viper_scores %>%
          dplyr::select(-std) %>%   
          spread(tf, avg) %>%
          data.frame(row.names = 1, check.names = FALSE) 

        summarized_viper_scores_df = summarized_viper_scores_df[ levels(Idents(anndataO_doro)), ]

        palette_length = 100
        my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

        
        message('saving csv as:', gsub('.pdf', '.csv', outfile))
        write.csv(summarized_viper_scores_df, file = gsub('.pdf', '.csv', outfile))
}



compute_TFactivities = function(DE_signatures, statistic='avg_logFC', pval = 'p_val_adj', cluster='cluster', gene='gene', statistic_th = 0){
    
    require(viper)
    require(reshape2)
    
    # Format input
    DE_signatures$statistic = DE_signatures[, statistic]
    DE_signatures$pval = DE_signatures[, pval]
    DE_signatures$cluster = DE_signatures[, cluster]
    DE_signatures$gene = DE_signatures[, gene]
    
    # Load regulons
    viper_gset = get(load('~/farm/gsea/genesets/dorotheav2-top10scoring_VentoLab20201111.rdata'))
    message('Analysizing ', length(viper_gset), ' TFs')
    
    # Compute activities
    TFact_holder = list()
    TFagreement_holder = list()
    for (cl_name in unique(DE_signatures$cluster)){
      message('\n', cl_name)

      cl_DEsig = subset(DE_signatures, cluster == cl_name)
      # Exclude probes with unknown or duplicated gene symbol
      cl_DEsig = subset(cl_DEsig, gene != "" )
      cl_DEsig = subset(cl_DEsig, ! duplicated(gene))
      rownames(cl_DEsig) = cl_DEsig$gene
        
      # Estimate z-score values for the GES. Cheeck VIPER manual for details
      mystatistic = matrix(cl_DEsig$statistic, dimnames = list(cl_DEsig$gene, 'statistic') )
      myPvalue = matrix(cl_DEsig$pval, dimnames = list(cl_DEsig$gene, 'pval') )
      mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(mystatistic))[, 1]
      mySignature = mySignature[order(mySignature, decreasing = T)]
      # Estimate TF activities
#       mrs = msviper(ges = mySignature, regulon = viper_gset, minsize = 4, ges.filter = F)
      mrs = msviper(ges = mystatistic[,1][order(mystatistic[,1], decreasing = T)], regulon = viper_gset, minsize = 3, ges.filter = F)
      cl_TFact = data.frame(Regulon = names(mrs$es$nes),
                                 cluster = cl_name,
                                 Size = mrs$es$size[ names(mrs$es$nes) ], 
                                 NES = mrs$es$nes, 
                                 p.value = mrs$es$p.value, 
                                 FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
      cl_TFact = subset(cl_TFact, Size < 200)
      cl_TFact = cl_TFact[ order(cl_TFact$p.value), ]
        
      cl_TFact$Regulon = as.character(cl_TFact$Regulon)
      cl_TFact$TF = sapply(strsplit(cl_TFact$Regulon, ' - '), head, 1)
      cl_TFact$TF = sapply(strsplit(cl_TFact$TF, '_'), head, 1)
        
      if( nrow(cl_TFact) > 0 )
        TFact_holder[[cl_name]] = cl_TFact
      
        
     if(nrow(cl_TFact) == 0)
        next()
       
      # Find TFs that are DEG and Dactive
      sharedTFs = intersect(subset(cl_DEsig, pval < 0.1 )$gene, subset(cl_TFact, FDR < 0.1)$TF)

      if( length(sharedTFs) > 0 ) {
        # Add TFs agreement info
        cl_TFagreement = data.frame(cluster = cl_name, TF = sharedTFs, stringsAsFactors = F)
        cl_TFagreement$TF_expression = 'upregulated'
        cl_TFagreement$TF_expression_statistic = cl_DEsig[sharedTFs, ]$statistic
        cl_TFagreement$TF_expression[ cl_TFagreement$TF_expression_statistic < statistic_th  ] = 'downregulated'
        cl_TFagreement$TF_expression_pval = cl_DEsig[sharedTFs, ]$pval
        cl_TFagreement$TF_activity = 'active_regulon'
        cl_TFagreement$TF_activity_FDR = cl_TFact$FDR[ match(cl_TFagreement$TF, cl_TFact$TF) ]
        cl_TFagreement$TF_activity_score = cl_TFact$NES[ match(cl_TFagreement$TF, cl_TFact$TF) ]
        cl_TFagreement$TF_activity[ cl_TFagreement$TF_activity_score < 0 ] = 'inactive_regulon'
        TFagreement_holder[[cl_name]] = cl_TFagreement
      }
        
    }

    # Merge TF activities
    TFactivities = melt(TFact_holder, id.vars = names(TFact_holder[[1]]))
    TFactivities = TFactivities[ order(TFactivities$NES, decreasing = T), -ncol(TFactivities) ]
    
    # Merge TF agreement
    TFagreement_holder = TFagreement_holder[ sapply(TFagreement_holder, length) > 0 ]
    if( length(TFagreement_holder)>0 ){
      TFagreement = melt(TFagreement_holder, id.vars = names(TFagreement_holder[[1]]))
      TFagreement$TF_expression_pval = signif(TFagreement$TF_expression_pval, 3)
      TFagreement$TF_activity_FDR = signif(TFagreement$TF_activity_FDR, 3)
      TFagreement = TFagreement[order(TFagreement$TF_activity_score, decreasing = T), ]
      TFagreement = TFagreement[order(TFagreement$cluster), ]
    }
    
    return(list(TFactivities = TFactivities, TFagreement = TFagreement))
    
}


plotHeatmap_TFact = function(TFact, TFs_of_interest, pdf_file = NA){
    df = subset(TFact$TFactivities, TF %in% TFs_of_interest)
    # build matrix of TF activity scores to plot
    df$value = df$NES
    NES = acast(df, TF~cluster, fill = 0)
    NES = NES[TFs_of_interest, ]
    labels = NES
    labels[] = ''
    for( i in rownames(labels) )
        for( j in colnames(labels) ){
            ag = subset(TFact$TFagreement, cluster == j & TF == i)
            if( nrow(ag) > 0 )
                labels[i , j] = '*'
            }
        
    pheatmap(t(NES), cellheight = 10, cellwidth = 10, cluster_rows = F, cluster_cols = F, display_numbers = t(labels),
             color = colorRampPalette(c(brewer.pal(n = 5, name = 'Blues')[4], "white", brewer.pal(n = 5, name = 'Reds')[4]))(50),
            filename = pdf_file)
}


plotHeatmap_TFexp = function(TFact, TFs_of_interest, pdf_file = NA){
    df = subset(TFact$TFagreement, TF %in% TFs_of_interest)
    # build matrix of TF activity scores to plot
    df$value = df$TF_expression_statistic
    FoldC = acast(df, TF~cluster, fill = 0)
    FoldC = FoldC[TFs_of_interest, ]
    labels = FoldC
    labels[] = ''
    for( i in rownames(labels) )
        for( j in colnames(labels) ){
            ag = subset(TFact$TFagreement, cluster == j & TF == i)
            if( nrow(ag) > 0 )
                labels[i , j] = '*'
            }
        
    pheatmap(t(FoldC), cellheight = 10, cellwidth = 10, cluster_rows = F, cluster_cols = F, display_numbers = t(labels),
             color = colorRampPalette(c(brewer.pal(n = 5, name = 'Blues')[4], "white", brewer.pal(n = 5, name = 'Reds')[4]))(50),
            filename = pdf_file)
}



compute_PathwayEnrichment = function(DE_signatures, statistic='avg_logFC', pval = 'p_val_adj', cluster='cluster', gene='gene', statistic_th = 0){
    
    require(viper)
    require(reshape2)
    
    # Format input
    DE_signatures$statistic = DE_signatures[, statistic]
    DE_signatures$pval = DE_signatures[, pval]
    DE_signatures$cluster = DE_signatures[, cluster]
    DE_signatures$gene = DE_signatures[, gene]
    
    # Load Pathways
#     viper_gset = get(load('~/farm/gsea/genesets/ovarianDevelopment_viper.rdata'))
#     viper_gset = append(x = viper_gset,  get(load('~/farm/gsea/genesets/stevant2019_viper.rdata')))
    
    viper_gset = get(load('~/farm/gsea/genesets/GObpUniprot_viper.rdata'))
    viper_gset = append(x = viper_gset,  get(load('~/farm/gsea/genesets/UniprotKeywords_viper.rdata')))
    viper_gset = append(x = viper_gset,  get(load('~/farm/gsea/genesets/ovarianDevelopment_viper.rdata')))
#     viper_gset = get(load('~/farm/gsea/genesets/UniprotKeywords_viper.rdata'))
#     viper_gset = viper_gset[ grep('signal', names(viper_gset), ignore.case = T) ]
#     viper_gset = viper_gset[ c(grep('KEGG', names(viper_gset), ignore.case = T),
#                               grep('REACTOME', names(viper_gset), ignore.case = T)) ]
    
#     viper_gset = get(load('~/farm/gsea/genesets/UniprotKeywords_viper.rdata'))
    message('Analysizing ', length(viper_gset), ' pathways')
    
    # Compute activities
    PathwayEn_holder = list()
    for (cl_name in unique(DE_signatures$cluster)){
      message('\n', cl_name)

      cl_DEsig = subset(DE_signatures, cluster == cl_name)
      # Exclude probes with unknown or duplicated gene symbol
      cl_DEsig = subset(cl_DEsig, gene != "" )
      cl_DEsig = subset(cl_DEsig, ! duplicated(gene))
      rownames(cl_DEsig) = cl_DEsig$gene
      if(nrow(cl_DEsig) < 5)
          next()
          
      # Estimate z-score values for the GES. Cheeck VIPER manual for details
      mystatistic = matrix(cl_DEsig$statistic, dimnames = list(cl_DEsig$gene, 'statistic') )
      myPvalue = matrix(cl_DEsig$pval, dimnames = list(cl_DEsig$gene, 'pval') )
      mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(mystatistic))[, 1]
      mySignature = mySignature[order(mySignature, decreasing = T)]
      # Estimate TF activities
#       mrs = msviper(ges = mySignature, regulon = viper_gset, minsize = 4, ges.filter = F)
      mrs = msviper(ges = mystatistic[,1][order(mystatistic[,1], decreasing = T)], regulon = viper_gset, minsize = 4, ges.filter = F)
      cl_PathwayEn = data.frame(Pathway = names(mrs$es$nes),
                                 cluster = cl_name,
                                 Size = mrs$es$size[ names(mrs$es$nes) ], 
                                 NES = mrs$es$nes, 
                                 p.value = mrs$es$p.value, 
                                 FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
      cl_PathwayEn = subset(cl_PathwayEn, Size < 200)
      cl_PathwayEn = cl_PathwayEn[ order(cl_PathwayEn$p.value), ]
        
      cl_PathwayEn$Pathway = as.character(cl_PathwayEn$Pathway)
        
      if( nrow(cl_PathwayEn) > 0 )
        PathwayEn_holder[[cl_name]] = cl_PathwayEn

    }

    # Merge TF activities
    PathwayEn = melt(PathwayEn_holder, id.vars = names(PathwayEn_holder[[1]]))
    PathwayEn = PathwayEn[ order(PathwayEn$NES, decreasing = T), -ncol(PathwayEn) ]
    
  
    return(PathwayEn)
    
}