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

make_traj_from_index = function(frames) {a=0*(1:20000); a[unlist(frames)]=1; a}

ab_analyse = function(filename, prefix){
  ab_data = read.table(filename, col.names = c('frame','ab.res'))
  ab_traj = ab_data %>% group_by(ab.res) %>% 
    summarise(frameIDs = list(frame)) %>% 
    mutate(traj = map(frameIDs, make_traj_from_index)) %>% 
    select(-frameIDs) %>%
    mutate(tframes = map_dbl(traj,~sum(unlist(.)))) %>%
    mutate(full_mean = tframes/20000) %>%
    filter(full_mean > .11) %>%
    left_join(ab_names)  %>%
    mutate(postsumm = map(traj, get_postsummary))  %>%   
    unnest(postsumm) %>% select(-traj)
  ab_traj %>% add_column(prefix=prefix)
}



int_analyse = function(filename, prefix){
  int_data = read.table(filename, col.names = c('frame','int.res'))
  int_traj = int_data %>% group_by(int.res) %>% 
    summarise(frameIDs = list(frame)) %>% 
    mutate(traj = map(frameIDs, make_traj_from_index)) %>% 
    select(-frameIDs) %>%
    mutate(tframes = map_dbl(traj,~sum(unlist(.)))) %>%
    mutate(full_mean = tframes/20000) %>%
    filter(full_mean > .11) %>%
    left_join(int_names) %>% 
    mutate(postsumm = map(traj, get_postsummary))  %>%   
    unnest(postsumm) %>% select(-traj)
  int_traj %>% add_column(prefix=prefix)
}


int_whole = int_analyse("contact_frame/trajectory/integrin.contacts.frame.txt",prefix='whole')
int_sc = int_analyse("contact_frame/trajectory/nobb.integrin.contacts.frame.txt",prefix='sc')


ab_whole = ab_analyse("contact_frame/trajectory/antibody.contacts.frame.txt",prefix='whole')
ab_sc = ab_analyse("contact_frame/trajectory/nobb.antibody.contacts.frame.txt",prefix='sc')



### ===================== END read trajectories, calc estimates ====================


rbind(ab_whole, ab_sc) %>% 
  ggplot(aes(x=ab.label, y=p1_50,color=prefix)) + 
  geom_point() + 
  facet_grid(.~ab.chain,scale='free_x') +
  geom_errorbar(aes(ymin=p1_05, ymax=p1_95), width=.3) +
  theme(axis.text.x=element_text(angle=90, hjust=1))  

rbind(int_whole, int_sc) %>% 
  ggplot(aes(x=int.label, y=p1_50,color=prefix)) + 
  geom_point() + 
  facet_grid(.~int.chain,scale='free_x',space='free') +
  geom_errorbar(aes(ymin=p1_05, ymax=p1_95), width=.3) +
  theme(axis.text.x=element_text(angle=90, hjust=1))  

rrr= rbind(int_whole, int_sc)

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
abcryo_sc = ab_read_cryo("contact_frame/cryo-em/nobb.cryo-em.antibody.contacts.txt",'cryo_SC')

intcryo_whole = int_read_cryo("contact_frame/cryo-em/cryo-em.integrin.contacts.txt",'cryo_whole')
intcryo_sc = int_read_cryo("contact_frame/cryo-em/nobb.cryo-em.integrin.contacts.txt",'cryo_SC')


redu = function(data) {
  data %>% 
    select(starts_with("p1"), 
           contains("label"), 
           contains("chain"), prefix)
  }


ab_tot = rbind(ab_whole %>% redu, ab_sc %>% redu) 

ggplot() + 
  geom_col(aes(x=ab.label, y=p1_50), fill='gray60', alpha=.2, position='dodge', 
           data=abcryo_whole %>% redu) +
  geom_point(
    aes(x=ab.label, y=p1_50, color=prefix, group=prefix),
    position=position_dodge(width = .5),size=2 , data=ab_tot) + 
  facet_grid(.~ab.chain,scale='free_x',space='free_x') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  geom_errorbar(aes(x=ab.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                width=.3,position=position_dodge(width = .5), data=ab_tot) +
  scale_color_brewer(palette = "Paired") 
  
  
  
  
int_tot = rbind(int_whole %>% redu, int_sc %>% redu) 
  
ggplot() + 
    geom_col(aes(x=int.label, y=p1_50), fill='gray60', alpha=.2, position='dodge', 
             data=intcryo_whole %>% redu) +
    geom_point(
      aes(x=int.label, y=p1_50, color=prefix, group=prefix),
      position=position_dodge(width = .5),size=2 , data=int_tot) + 
    facet_grid(.~int.chain,scale='free_x', space='free_x') +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    geom_errorbar(aes(x=int.label, ymin=p1_05, ymax=p1_95,color=prefix, group=prefix), 
                  width=.3,position=position_dodge(width = .5), data=int_tot) +
    scale_color_brewer(palette = "Paired") 
  
