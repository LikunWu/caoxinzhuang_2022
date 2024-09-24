library(tidyverse)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

bac_bins_tree <- read.tree('bins_tree.nwk')
bac_bin_asv <- read.csv('bacteria.Bin_TopTaxon.csv',header = T)
bac_binstoprocess <- read.csv('bacteria.BinContributeToProcess_EachGroup.csv',header = T)
bac_rt_process <- read.csv('bacteriaiCAMP_outProcessImportance_EachGroup.csv',header = T)
bac_bin_phylum1 <- read.csv('bacteria.Taxon_Bin.csv',header = T)

col_icamp <- c('#2A9E89','#FFCB4E','#DA546D','#F7B7AE','#E56F51')
CAOXINZHUANG2 <- c('#B44847','#00A1D6')
sunsx <- c('#A8DDE3','#F7E1AF','#EEBCBF','#466E88','#DCB46C','#D6585B')
#bacteria communities assembly process under H and R

bac_icamp_total_data <- bac_rt_process %>%
  filter(Group != 'Remove_vs_Turnover') %>%
  pivot_longer(!c(1:3),names_to = 'process',values_to = 'value') %>%
    group_by(Group) %>%
    mutate(ymax = cumsum(value)) %>%
    group_by(Group) %>%
    mutate(ymin = c(0,head(ymax,n = -1)))%>%
  mutate(labelposition = (ymax+ymin)/2)

####Figure 4A####
bac_icamp_total_p2 <- ggplot(bac_icamp_total_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=process)) +
  geom_rect() +
  scale_fill_manual( 
    limits=c("HoS","HeS","HD","DR","DL"),
    values=col_icamp,
    guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2,
                       override.aes=list(alpha=0.8)))+
  geom_text( x=3, aes(y=labelposition, label=paste(round(value*100,digits = 2),'%'), 
                      color=process), size=6) +
  coord_polar(theta="y") +
  scale_color_manual(limits=c("HoS","HeS","HD","DR","DL"),
                     values=col_icamp,
                     guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2,
                                        override.aes=list(alpha=0.8))) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  facet_grid(.~Group)
##########################################################################
#####################Figure 4A finished######################################
##########################################################################

#####bacteria bins #####
# extract sign for each node 
# delete quote sign in the tree tip label
bac_bins_tree$tip.label <- str_replace_all(bac_bins_tree$tip.label,'[^ASV_\\d]','')
# annotate phylum to bins
# caution! some bins are composed ASVs from more than one phylum
bacteria_bin_unique <-  bac_bin_phylum1 %>% 
  select(Bin,Phylum) %>%
  filter(Phylum!='') %>%
  group_by(Bin) %>%
  unique()%>%
  mutate(phylum2 = ifelse(isUnique(Bin),Phylum,'Uncer_Bac_Phy')) %>%
  left_join(bac_bin_asv,by = 'Bin') 
bacteria_bin_unique2 <- select(bacteria_bin_unique,c(phylum2,TopTaxonID)) %>%
  as.data.frame()

# annotate phylum to tree branch
bac_bins_phylum_group <- list(Acidobacteriota = filter(bacteria_bin_unique2,
                                                        phylum2 == 'Acidobacteriota') %>% select(TopTaxonID) %>% as.list()%>% .[[1]],
                              Actinobacteriota = filter(bacteria_bin_unique2,
                                                     phylum2 == 'Actinobacteriota') %>% select(TopTaxonID)%>% as.list() %>% .[[1]],
                              Chloroflexi  = filter(bacteria_bin_unique2,
                                                         phylum2 == 'Chloroflexi') %>% select(TopTaxonID)%>% as.list()%>% .[[1]],
                              Proteobacteria  = filter(bacteria_bin_unique2,
                                                        phylum2 == 'Proteobacteria') %>% select(TopTaxonID)%>% as.list()%>% .[[1]],
                              Unspecified_Phylum = filter(bacteria_bin_unique2,
                                              phylum2 == 'Uncer_Bac_Phy') %>% select(TopTaxonID)%>% unique()%>%as.list()%>% .[[1]],
                              Others  = filter(bacteria_bin_unique2,
                                                       phylum2 %in% c('Armatimonadota','Bdellovibrionota',
                                                                      'Cyanobacteria','Firmicutes',
                                                                      'Gemmatimonadota','Nitrospirota','Patescibacteria',
                                                                      'Planctomycetota','SAR324 clade(Marine group B)','Verrucomicrobiota')) %>%
                                select(TopTaxonID)%>% as.list()%>% .[[1]])
bac_bins_group <- groupOTU(bac_bins_tree,bac_bins_phylum_group) # add phyla data to tree branch

# prepare phyla for other data
bacteria_bin_unique3 <- bacteria_bin_unique2  %>%
  mutate(phylum3 = ifelse(phylum2 %in% c('Armatimonadota','Bdellovibrionota',
                                         'Cyanobacteria','Firmicutes',
                                         'Gemmatimonadota','Nitrospirota','Patescibacteria',
                                         'Planctomycetota','SAR324 clade(Marine group B)','Verrucomicrobiota'),
                          'Others',phylum2))  %>%#annotate to color rings 2-4
  select(Bin,TopTaxonID,phylum3) %>%
  unique()

bacteria_bin_unique3$phylum3 %>% list() %>% table()
# ecological processes to every phylum
phylum_process_bac_relative <- bac_binstoprocess[,-c(1:2)] %>%
  filter(Group != 'Remove_vs_Turnover') %>%
  pivot_longer(!c(Group,Process),names_to = 'bin',values_to = 'value') %>%
  na.omit() %>%
  mutate(Bin = str_replace_all(bin, 'bin' , 'Bin')) %>%
  left_join(bacteria_bin_unique3,by = 'Bin') %>%
  group_by(phylum3,Process) %>%
  summarise(sumvalue = sum(value)) %>%
  group_by(phylum3) %>%
  summarise(relative_process = sumvalue/sum(sumvalue)) %>%
  mutate(process = rep(c('DL','DR','HD','HeS','HoS')))

# ecological processes to each bin
bac_bin_HoSDR_relative <- bac_binstoprocess[,-c(1:2)] %>%
  filter(Group != 'Remove_vs_Turnover') %>%
  pivot_longer(!c(Group,Process),names_to = 'bin',values_to = 'value') %>%
  na.omit() %>%
  filter(Process %in% c('HoS','DR')) %>%
  pivot_wider(names_from = 'Group',values_from = 'value') %>%
  group_by(Process,bin) %>%
  mutate(relative = Turnover - Remove) %>%
  mutate(Bin = str_replace_all(bin , 'bin' , 'Bin')) %>% # HoS和DR受到秸秆添加的相对影响
  left_join(bac_bin_asv,by = 'Bin')  

####relative abundance of each bins####
bac_bin_relative_abundance <- bac_bin_phylum1 %>%
  group_by(Bin) %>%
  summarise(abundance = sum(TaxonRelativeAbundance)) %>%
  arrange(desc(abundance)) %>%
  left_join(bacteria_bin_unique3,by = 'Bin') 

bac_bin_phyla_relative_abundance <- bac_bin_relative_abundance %>%
  group_by(phylum3) %>%
  summarise(relativeabundance = sum(abundance))

bac_binstoprocess_long <- bac_binstoprocess[,-c(1:2)] %>%
  filter(Group=='Remove_vs_Turnover') %>%
  pivot_longer(!c(Group,Process),names_to = 'bin',values_to = 'value') %>%
  na.omit() %>%
  group_by(bin) %>%
  summarise(relative_process = value / sum(value)) %>% #转换为每个过程的相对值
  mutate(process = rep(c('HeS','HoS','DL','HD','DR'),n = length(bac_bins_tree$tip.label))) %>%
  mutate(Bin = str_replace_all(bin , 'bin' , 'Bin')) 

bac_binstoprocess_long %>%
  filter(relative_process>=0.5) %>%
  group_by(process) %>%
  summarise(count = table(bin)) %>%
  dplyr::count(process) %>%
  mutate(proportion = n/179)

bac_bins_1st_cir <- data.frame(ID = bac_bins_tree$tip.label) %>%
  left_join(bac_bin_asv[,c(1,3)], by =c ('ID'='TopTaxonID')) %>%
  left_join(bac_binstoprocess_long,by = 'Bin') #relative importance for each process in each bin
bac_bins_2nd_cir <-data.frame(ID = bac_bins_tree$tip.label) %>% 
  left_join(bac_bin_asv[,c(1,3)], by =c ('ID'='TopTaxonID')) %>%
  left_join(bac_bin_relative_abundance,by = 'Bin')
bac_bins_3rd_cir_HoS <- data.frame(ID = bac_bins_tree$tip.label) %>%
  left_join(select(bac_bin_HoSDR_relative,c(Process,relative,TopTaxonID,Bin)), by =c ('ID'='TopTaxonID')) %>%
  left_join(bacteria_bin_unique3,by = 'Bin') %>%
  filter(Process == 'HoS')
bac_bins_4th_cir_DR <- data.frame(ID = bac_bins_tree$tip.label) %>%
  left_join(select(bac_bin_HoSDR_relative,c(Process,relative,TopTaxonID,Bin)), by =c ('ID'='TopTaxonID')) %>%
  left_join(bacteria_bin_unique3,by = 'Bin') %>%
  filter(Process == 'DR')

bac_bin_icamp <-
  ggtree(bac_bins_group,aes(col = group),
       layout="circular", branch.length='none', size=0.2) + 
  geom_treescale(x=-15, y=0, fontsize=0.1, linesize=0.1)+
    # geom_tiplab(size=3)+
  scale_color_manual( values = sunsx)+
  new_scale_fill()+
  geom_fruit(data = bac_bins_1st_cir,
             geom = geom_bar,
             mapping = aes(x = relative_process,y = ID,fill= process),alpha = 0.7,
             stat = 'identity',
             orientation = 'y', width = 1,position="auto",
             pwidth=1) +
  scale_fill_manual( 
    limits=c("HoS","HeS","HD","DR","DL"),
    values=col_icamp,
    guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2,
                       override.aes=list(alpha=0.8)))+
  new_scale_fill()+
  geom_fruit(data = bac_bins_2nd_cir,
             geom = geom_bar,
             mapping = aes(x = log10(100*abundance+1),y = ID,fill= phylum3),alpha = 0.7,
             stat = 'identity',
             orientation = 'y', width = 1,position="auto",offset = 0.3,
             pwidth=0.3) +
  scale_fill_manual(values = sunsx)+
  new_scale_fill()+
  geom_fruit(data = bac_bins_3rd_cir_HoS,
             geom = geom_col,
             mapping = aes(x = relative ,y = ID,fill= phylum3),
             stat = 'identity',
             orientation = 'y', width = 0.8,position="auto",offset = 0.3,
             pwidth=0.3)+
  scale_fill_manual(values = sunsx)+
  new_scale_fill()+
  geom_fruit(data = bac_bins_4th_cir_DR,
             geom = geom_bar,
             mapping = aes(x = relative,y = ID,fill= phylum3),
             stat = 'identity',
             orientation = 'y', width = 0.8,position="auto",offset = 0.3,
             pwidth=0.3)+
  scale_fill_manual(values = sunsx)
bac_bin_icamp

##################################################################
##########Figure 5A finished######################################
##################################################################

####ecological processes contribute to each phylum####
 phylum_process_bac_relative %>%
  ggplot()+
  geom_col(aes(x = phylum3,y = relative_process ,fill = process))+
  scale_fill_manual(
    limits=c("HoS","HeS","HD","DR","DL"),
    values=col_icamp,
    guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2,
                       override.aes=list(alpha=0.8))
  )+
  labs(x = 'Phyla',y = 'Relative contribution')+
  theme_NMDS+
  theme(legend.position = '')

bac_bin_phyla_relative_abundance_p <- bac_bin_phyla_relative_abundance %>%
  mutate(aaa = 1) %>%
  ggplot()+
  geom_bar(aes(x = aaa,y = relativeabundance ,fill = phylum3  ),stat = 'identity',
           position = 'stack')+
  labs(x = NULL,y = 'Relative abundance')+
  theme_NMDS+
  scale_fill_manual(values = sunsx)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
bac_phyla_process_p <- ggarrange(phylum_process_bac_relative_p,bac_bin_phyla_relative_abundance_p,widths = c(8,2),heights = 4)

#### Phylum contribution to HoS and DR####
bac_phylum_contribution_to_process <- bac_binstoprocess %>%
  filter(Group == 'Remove_vs_Turnover') %>%
  select(!c(Method,GroupBasedOn,Group)) %>%
  pivot_longer(!Process,names_to = 'bin',values_to = 'value') %>%
  group_by(Process) %>%
  mutate(relative = value/sum(value)) %>%
  mutate(Bin = str_replace_all(bin, 'bin','Bin')) %>%
  left_join(bacteria_bin_unique3,by = 'Bin') %>%
  group_by(Process,phylum3) %>%
  summarise(phylum_relative = sum(relative))%>%
  filter(Process %in% c('HoS','DR')) #importance of different processes contribute to each bins

bac_phylum_contribution_to_process
bac_phylums_contribution_p <- ggplot(bac_phylum_contribution_to_process) +
   geom_col(aes(x = Process,y = phylum_relative,fill = phylum3))+
  scale_fill_manual(values = sunsx)+
  labs(y = 'Phyla contribution')+
  theme_classic()+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 23),
        axis.title.x = element_blank())
###########################Figure S7A###########################