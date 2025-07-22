###Author: Shihan Li
comb <- readRDS('231201_zheng_final.rds')

########### Fig.3e/f Umaps of Zheng dataset with TEX and TRM signature overlaid ###########

FeaturePlot(comb, features =c('new_trm_sig_1', 'new_tex_sig_1'), 
            split.by = 'isCancer', min.cutoff = 0, max.cutoff = 'q99', order=T) & 
    coord_fixed(ratio=0.8) & 
    scale_colour_continuous(plotto_color_gradient_blue_red)

FeaturePlot(comb, features =c('new_trm_sig_1'), order=T, min.cutoff = 0, max.cutoff = 'q99') + 
    
    for(feature in c('new_trm_sig_1', 'new_tex_sig_1')) {
        
        p1 <- FeaturePlot(comb[,comb$isCancer == 'Healthy'], features =feature, min.cutoff = 0, max.cutoff = 'q99', order=T) + 
            coord_fixed(ratio=1) + 
            scale_colour_viridis_c(option = 'E') + ggtitle('Healthy')
        
        p2 <- FeaturePlot(comb[,comb$isCancer == 'Tumour'], features =feature, min.cutoff = 0, max.cutoff = 'q99', order=T) + 
            coord_fixed(ratio=1) + 
            scale_colour_viridis_c(option = 'E') + ggtitle('Tumour')
        
        p3 <- FeaturePlot(comb[,comb$isCancer == 'Healthy'], features =feature, min.cutoff = 0, max.cutoff = 'q99', order=T) + 
            coord_fixed(ratio=1) + scale_colour_gradientn(colours = plotto_color_gradient_blue_red)+ ggtitle('Healthy')
        
        p4 <- FeaturePlot(comb[,comb$isCancer == 'Tumour'], features =feature, min.cutoff = 0, max.cutoff = 'q99', order=T) + 
            coord_fixed(ratio=1) + scale_colour_gradientn(colours = plotto_color_gradient_blue_red) + ggtitle('Tumour')
        
        p5 <- (p1|p2)/(p3|p4)
        ggsave(p5, filename = paste0('plots/231123_', feature,'.pdf'), dpi = 300, scale = 5)
        
    }


######### Figure 3g, Bar plots with relative % of TEX, TEM, TRM from each tumour type ########

p1 <- ggplot(comb@meta.data,aes(x=cancerType, fill=metaclusters_coarse)) + geom_bar(position='fill', colour='black') + 
    facet_grid(~isCancer) & theme(legend.position = 'none',
                                  panel.background = element_blank(),
                                  panel.border = element_rect(colour='black', fill=NA),
                                  axis.line = element_blank())

p1 <- ggplot(comb@meta.data,aes(x=cancerType, fill=metaclusters_coarse)) + geom_bar(position='fill', colour='black') & theme(legend.position = 'none',
                                                                                                                             panel.background = element_blank(),
                                                                                                                             panel.border = element_rect(colour='black', fill=NA),
                                                                                                                             axis.line = element_blank())
p1
ggsave(p1, filename = paste0('plots/231123_proportions_byCelltype.pdf'), dpi = 300, scale = 3)



######### Figure 4a, heatmaps TCR clones ########


temp <- comb@meta.data %>% group_by(CDR3.A1, CDR3.B1) %>% mutate(publicity = n_distinct(patient),
                                                                 functional_clone_id = dplyr::cur_group_id(),
                                                                 clone_count = n()) %>% as.data.frame()
row.names(temp) <- temp$cellID
comb@meta.data <- temp
1

meta <- read.csv('tcr_analysis/240530_full_dataset_burn.csv', na.strings = '')

meta <- meta %>% group_by(cdr3_a_nucseq, cdr3_b_nucseq, .drop = T) %>% mutate(real_clone_id = cur_group_id(),
                                                                              real_clone_count = n())
meta[meta$v_a_gene == 'TRAV1-2*01' & meta$j_a_gene %in% c('TRAJ12*01', 'TRAJ20*01', 'TRAJ33*01'),'celltype'] <- 'MAIT'
meta[meta$celltype %in% c('TCM', 'NAIVE'),'celltype'] <- 'NAIVE_TCM'
meta[meta$celltype %in% c('gdT', 'ISG_CD8', 'TPEX'),'celltype'] <- 'other'

meta[meta$origin == 'burn' & meta$cancerType == 'HC','origin'] = 'burn_HC'
meta[meta$origin == 'burn' & meta$cancerType == 'BC','origin'] = 'burn_BC'


get_order <- function(df, celltype_order, group_id){
    to_return = c()
    for(celltype in celltype_order){
        df = df[!df[[group_id]] %in% to_return,]
        temp = df[df[[group_id]] == celltype,]
        temp = arrange(temp, prop)
        temp = temp[temp$prop > 0.5,]
        to_return = append(to_return, unique(temp[['real_clone_id']]))
    }
    return(to_return)
}

for_loop = c('zheng', 'burn_BC', 'burn_HC')
for(var in for_loop){
    clones <- meta[meta$origin == var,] %>% arrange(desc(real_clone_count))
    
    
    clone <- unique(clones$real_clone_id)[1:100]
    clones <- clones[clones$real_clone_id %in% clone,]
    clones$celltype <- factor(clones$celltype)
    
    
    temp <- clones %>% group_by(real_clone_id,celltype,.drop = FALSE) %>%
        summarise(count = n()) %>%
        mutate(prop = count / sum(count)) 
    
    
    z <- get_order(temp, celltype_order = c("TRM", "TEX", "TEM", "TEMRA", "NAIVE_TCM","MAIT","other"), 'celltype')
    
    p1 <- ggplot(temp, 
                 aes(x=factor(real_clone_id, levels = z), y=factor(celltype, 
                                                                   levels = c("TRM", "TEX", "TEM", "TEMRA", "NAIVE_TCM","MAIT","other")), 
                     fill=prop)) + 
        geom_tile() + scale_fill_gradientn(colours = plotto_color_gradient_blue_red) + ylab('Celltype') + xlab('Top 100 expanded clones') + 
        scale_y_discrete(expand = c(0,0)) + 
        ggtitle(var) +
        theme(axis.text.x = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour='black', fill=NA),
              axis.line = element_blank())
    ggsave(p1, filename = paste0('plots/240618_heatmap_100_', var,'.pdf'), dpi=300, scale=3)
}

