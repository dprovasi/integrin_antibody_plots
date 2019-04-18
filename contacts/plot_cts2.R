setwd("~/Desktop/integrin_antibody_figures/contacts")


### ===================== read names, recode chain name ====================

read_residues = function(filename, prefix=NULL, chain_recode=NULL){
  ### implement recode for chain names!
  if(is.list(chain_recode) ){
    recode_chain = function(data) {
      data %>% mutate(chain = recode(chain, !!!chain_recode))}
  } else { recode_chain = function(data) data }
  
  read.table(filename, header = TRUE) %>% 
    separate(resid, c("resname","resID"), sep =3) %>%
    mutate(res= aa321(resname))   %>% 
    recode_chain() %>%
    mutate(label = str_c(chain, str_c(res, resID),sep = "/")) %>%
    mutate(label = fct_reorder(label, order(order(chain,residue)))) %>%
    setNames(paste0(prefix, names(.))) 
}

ab_names = read_residues(
  "contact_frame/antibody.names.txt",'ab.', 
  chain_recode = list(C='H',D='L')) %>%
  mutate(ab.res1 = ab.res, ab.res=as.numeric(ab.residue))

int_names = read_residues(
  "contact_frame/integrin.names.txt",'int.',
  chain_recode = list(A='aIIb',B='b3')) %>%
  mutate(int.res1 = int.res, int.res=as.numeric(int.residue))
### ===================== END read names, recode chain name ====================




### ===================== read trajectories, calc estimates ====================
trj_from_sparse = function(frames) {a=0*(1:20000); a[unlist(frames)]=1; a}

ab_analyse = function(filename, prefix, make_traj_from_index= trj_from_sparse ){
  ab_data = read.table(filename, col.names = c('frame','ab.res'))
  ab_traj = ab_data %>% group_by(ab.res) %>% 
    summarise(frameIDs = list(frame)) %>% 
    mutate(traj = map(frameIDs, make_traj_from_index)) %>% 
    select(-frameIDs) %>%
    mutate(tframes = map_dbl(traj,~sum(unlist(.)))) %>%
    mutate(full_mean = map_dbl(traj,~mean(unlist(.)))) %>%
    mutate(postsumm = map(traj, get_postsummary))  %>%   
    unnest(postsumm) %>% select(-traj)
  ab_traj %>% add_column(prefix=prefix)
}


int_analyse = function(filename, prefix,make_traj_from_index= trj_from_sparse){
  int_data = read.table(filename, col.names = c('frame','int.res'))
  int_traj = int_data %>% group_by(int.res) %>% 
    summarise(frameIDs = list(frame)) %>% 
    mutate(traj = map(frameIDs, make_traj_from_index)) %>% 
    #select(-frameIDs) %>%
    mutate(tframes = map_dbl(traj,~sum(unlist(.)))) %>%
    mutate(full_mean = tframes/20000) %>%
    mutate(postsumm = map(traj, get_postsummary))  %>%   
    unnest(postsumm) %>% select(-traj)
  int_traj %>% add_column(prefix=prefix)
}


int_whole = int_analyse("contact_frame/trajectory/integrin.contacts.frame.txt",prefix="bb+sc")
int_sc = int_analyse("contact_frame/trajectory/nobb.integrin.contacts.frame.txt",prefix='sc')
ab_whole = ab_analyse("contact_frame/trajectory/antibody.contacts.frame.txt",prefix="bb+sc")
ab_sc = ab_analyse("contact_frame/trajectory/nobb.antibody.contacts.frame.txt",prefix='sc')

### ===================== END read trajectories, calc estimates ====================


### ===================== read cryo_EM data ====================

ab_read_cryo  = function(filename, prefix) {
  read.table(filename, col.names = c('frame','ab.res')) %>% 
    select(-frame) %>%
    left_join(ab_names) %>% 
    mutate(p1_50 = 1, p1_05=NA, p1_95=NA) %>%
    mutate(prefix=prefix)
}

int_read_cryo  = function(filename, prefix) {
  read.table(filename, col.names = c('frame','int.res')) %>% 
    select(-frame) %>%
    left_join(int_names) %>% 
    mutate(p1_50 = 1, p1_05=NA, p1_95=NA) %>%
    mutate(prefix=prefix)
}

abcryo_whole = ab_read_cryo("contact_frame/cryo-em/cryo-em.antibody.contacts.txt",'cryo_whole')
intcryo_whole = int_read_cryo("contact_frame/cryo-em/cryo-em.integrin.contacts.txt",'cryo_whole')


### filter out from the traj data the ones that are present in the cryo

redu = function(data) {
  data %>% 
    select(starts_with("p1"), 
           contains("label"), 
           contains("chain"), prefix)}

tt_int_whole = int_whole %>% 
  left_join(intcryo_whole %>% select(int.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(int_names)
tt_int_sc = int_sc %>% 
  left_join(intcryo_whole %>% select(int.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(int_names)

tt_ab_whole = ab_whole %>% 
  left_join(abcryo_whole %>% select(ab.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(ab_names)
tt_ab_sc = ab_sc %>% 
  left_join(abcryo_whole %>% select(ab.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(ab_names)

int_tot = rbind(tt_int_whole %>% redu, tt_int_sc %>% redu)
ab_tot = rbind(tt_ab_whole %>% redu, tt_ab_sc %>% redu) 

ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=ab.label, y=p1_50, color=prefix, group=prefix),
    position=position_dodge(width = .5),size=2 , data=ab_tot) +
  facet_grid(.~ab.chain,scale='free_x',space='free_x') +
  geom_errorbar(aes(x=ab.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5), data=ab_tot) +
  scale_color_brewer(palette = "Paired") +
  geom_col(aes(x=ab.label, y=p1_50), fill='gray60', alpha=.2, position='dodge', 
           data=abcryo_whole %>% redu) 
ggsave("contacts_ab.pdf", width = 9, height = 4)



ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=int.label, y=p1_50, color=prefix, group=prefix), 
    position=position_dodge(width = .5), size=2 , data=int_tot) + 
  facet_grid(.~int.chain,scale='free_x', space='free_x') +
  geom_errorbar(aes(x=int.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5), data=int_tot) +
  scale_color_brewer(palette = "Paired") +
  geom_col(aes(x=int.label, y=p1_50), fill='gray60', alpha=.2, position='dodge',
           data=intcryo_whole %>% redu)
ggsave("contacts_int.pdf", width = 8, height = 4)




  
####------ simple version

rbind(abcryo_whole %>% redu %>% mutate(prefix=recode(prefix, cryo_whole='cryo')), 
      ab_tot %>% filter(prefix == 'sc')) %>% ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=ab.label, y=p1_50, color=prefix, group=prefix),
    position=position_dodge(width = .5),size=2) +
  facet_grid(.~ab.chain,scale='free_x',space='free_x') +
  geom_errorbar(aes(x=ab.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5)) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "residue", y='probability')
ggsave("simple-contacts_ab.pdf", width = 9, height = 4)


rbind(intcryo_whole %>% redu %>% mutate(prefix=recode(prefix, cryo_whole='cryo')), 
      int_tot %>% filter(prefix == 'sc')) %>% ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=int.label, y=p1_50, color=prefix, group=prefix),
    position=position_dodge(width = .5),size=2) +
  facet_grid(.~int.chain,scale='free_x',space='free_x') +
  geom_errorbar(aes(x=int.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5)) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "residue", y='probability')
ggsave("simple-contacts_int.pdf", width = 8, height = 4)




###### -------------------------------- #####
### our clusters start from 1

f2c = read.table("cluster/c_num_vs_frame",col.names = c('frame','cluster')) %>%
  mutate(cluster = cluster + 1 )

## check populations
{
  cluster_pop = read.table("cluster/fraction.txt", col.names = 'population') %>% 
    add_column(cluster = (1:20))
  f2c %>% group_by(cluster) %>% 
    summarise(pop=n()) %>% 
    mutate(fraction=pop/10000) %>%
    left_join(cluster_pop) %>% 
    ggplot(aes(x=fraction, y=population)) + geom_point() + geom_abline(slope=1)
}

trj_from_cluster = function(clusters, mapf2c=f2c) {
  ### mapf2c is a table with frame and cluster columns
  a=0*(1:10000); ## we have trajs of 10k frames
  frames = unlist(mapf2c %>% filter(cluster %in% clusters) %>% select(frame))
  a[unlist(frames)]=1; 
  rep(a,each=2) ## to rematch 20k trajectories and keep the same lag time
  }


int_cluster = int_analyse(
  "cluster/integrin.cluster.contacts.frame.txt",
  prefix="cluster", 
  make_traj_from_index = trj_from_cluster)

#tmat2 = list(zz= 75, uu =1921, zu=1, uz=1)
#post_summary(tmat2, nbin=2000, tolerance=1e-2) 
#
#post_summary2(tmat2, nbin=1000, tolerance=1e-2) 

tt_int_cluster = int_cluster %>% 
  left_join(intcryo_whole %>% select(int.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(int_names)

reprefix = function(data){
  data %>% redu %>% mutate(prefix=recode(prefix, 'bb+sc'='traj:bb+sc', sc='traj:sc', cluster='clust')) %>%
    mutate(prefix=factor(prefix, levels=c('traj:bb+sc','traj:sc','clust')))}


int_tot2 = rbind(tt_int_whole %>% reprefix, tt_int_sc %>% reprefix, tt_int_cluster%>%reprefix) 

ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=int.label, y=p1_50, color=prefix, group=prefix), 
    position=position_dodge(width = .5), size=2 , data=int_tot2) + 
  facet_grid(.~int.chain,scale='free_x', space='free_x') +
  geom_errorbar(aes(x=int.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5), data=int_tot2) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values=c('#a6cee3', '#1f78b4', '#6a3d9a')) +
  geom_col(aes(x=int.label, y=p1_50), fill='gray60', alpha=.2, position='dodge',
           data=intcryo_whole %>% redu)
ggsave("contacts_int.pdf", width = 8, height = 4)



ab_cluster = ab_analyse(
  "cluster/antibody.cluster.contacts.frame.txt",
  prefix="cluster", 
  make_traj_from_index = trj_from_cluster)

tt_ab_cluster = ab_cluster %>% 
  left_join(abcryo_whole %>% select(ab.res, has_contact=p1_50)) %>%
  filter(has_contact==1 | full_mean>.11) %>% left_join(ab_names)

ab_tot2 = rbind(tt_ab_whole %>% reprefix, tt_ab_sc %>% reprefix, tt_ab_cluster %>% reprefix) 

ggplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_point(
    aes(x=ab.label, y=p1_50, color=prefix, group=prefix),
    position=position_dodge(width = .5),size=2 , data=ab_tot2) +
  facet_grid(.~ab.chain,scale='free_x',space='free_x') +
  geom_errorbar(aes(x=ab.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5), data=ab_tot2) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values=c('#a6cee3', '#1f78b4', '#6a3d9a')) +
  geom_col(aes(x=ab.label, y=p1_50), fill='gray60', alpha=.2, position='dodge', 
           data=abcryo_whole %>% redu) 
ggsave("contacts_ab.pdf", width = 9, height = 4)


#### ===== ===== ===== plot cc ========== ===== ===== 

int_data = read.table(
  "cluster/integrin.cluster.contacts.frame.txt", 
  col.names = c('cluster','int.res')) %>% left_join(int_names) 

ab_data = read.table(
  "cluster/antibody.cluster.contacts.frame.txt", 
  col.names = c('cluster','ab.res')) %>% left_join(ab_names) 

### order clusters with dendrogram -- 
order_clusters=function(data){
  mm = data %>% select(cluster, ends_with("label")) %>% drop_na() %>% 
      add_column(contact=1) %>% spread(cluster, contact) %>% select(-ends_with("label"))
  mm = as.matrix(mm)
  mm[is.na(mm)]=0
  ht=heatmap(mm, scale='none') 
  data %>% mutate(clusterID = factor(cluster, levels=ht$colInd))
}

##### drop_na since one of the residue names is missing from the files!!!

int_data %>% order_clusters() %>% drop_na() %>% 
  left_join(cluster_pop) %>% 
  ggplot(aes(y=clusterID, x=int.label, fill=population)) + 
  facet_grid(.~int.chain, scale='free_x',space='free_x') +
  geom_tile(colour="white",size=0.25) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_distiller(palette = "RdPu",trans = "reverse")
ggsave("int_cluster_description.pdf", width = 9, height = 4)


ab_data %>% order_clusters() %>% drop_na() %>%
  left_join(cluster_pop) %>% 
  ggplot(aes(y=clusterID, x=ab.label, fill=population)) + 
  facet_grid(.~ab.chain, scale='free_x',space='free_x') +
  geom_tile(colour="white",size=0.25) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_distiller(palette = "RdPu",trans = "reverse")
ggsave("ab_cluster_description.pdf", width = 9, height = 4)






