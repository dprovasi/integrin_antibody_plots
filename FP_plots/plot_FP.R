require(tidyverse)
require(bio3d) ## just for aa321
require(RColorBrewer)
setwd("~/Desktop/integrin_antibody_figures/FP_plots")


### --------------------------- prepare data ---------------------------------------
### read the residue names and change column names
### change to 1 letter repres

read_residues = function(filename, prefix=NULL){
  read.table(filename, header = TRUE) %>% 
    separate(resid, c("resname","resID"), sep =3) %>%
    mutate(res= aa321(resname))   %>% 
    mutate(label = str_c(chain, str_c(res, resID),sep = "/")) %>%
    setNames(paste0(prefix, names(.))) 
}

int_names = read_residues("integrin_resid.txt", prefix='int.')
ab_names = read_residues("antibody_resid.txt", prefix='ab.')


### read interaction names
fp_names = read.table("fingerprint_names.txt", 
                      header = FALSE, 
                      col.names = 'interaction') %>% 
  filter(interaction != "Wat_Hbond")

### prepare lines for dataset
n_int = nrow(int_names)
n_ab = nrow(ab_names)
n_fp = nrow(fp_names)
df_fpnames = fp_names[rep(seq_len(n_fp),n_int*n_ab), , drop=FALSE]
df_abnames = ab_names[rep(rep(seq_len(n_ab),each=n_fp),n_int), ]
df_intnames = int_names[rep(seq_len(n_int),each=n_fp*n_ab), ]

## read the fingerprint data for all frames

## calculate mean and sd for full set from frame matrix--
sift = read.table("sift_array.500ns",header=FALSE)
hascontact = (colSums(sift)>20)
sift_contact = sift[,hascontact] 
trajs=lapply(seq_along(which(hascontact)), function(i) sift_contact[,i])

data3 = cbind(
  df_fpnames[hascontact,, drop=FALSE],
  df_abnames[hascontact,],
  df_intnames[hascontact,]) %>% 
  add_column(traj=trajs)  %>% 
  add_column(full_mean = colMeans(sift_contact)) %>%
  rowid_to_column("ID")

  

#### ========= Bayesian Markov Model estimation ==========
tmat = function(traj,stride=1) {
  traj_stride = traj[seq(1, length(traj)-stride, stride)]
  traj1=head(traj_stride,-1)
  traj2=tail(traj_stride,-1)
  zz = sum(traj1==0 & traj2==0)
  zu = sum(traj1==0 & traj2==1)
  uz = sum(traj1==1 & traj2==0)
  uu = sum(traj1==1 & traj2==1)
  list(zz=zz,uu=uu,zu=zu,uz=uz)
}

ldp01 = function(p01,tm){(tm$zu-1)*log(p01)+(tm$zz-1)*log(1-p01) -lbeta(tm$zu,tm$zz)}
ldp10 = function(p10,tm){(tm$uz-1)*log(p10)+(tm$uu-1)*log(1-p10) -lbeta(tm$uz,tm$uu)}

post_summary =function(tm, 
                       cutoff=1e-10, 
                       verbose=FALSE, 
                       nbin=1000, tolerance=1e-2){
  ## calculates the summary in pi1 on the posterior distribution
  eps = 1e-8
  if(tm$uu==0 & tm$zu==0) {
    res = c(0,0,0)
    names(res) = c("p1_05","p1_50","p1_95")
    return(res)
  }
  if(tm$zz==0 & tm$uz==0) {
    res = c(1,1,1)
    names(res) = c("p1_05","p1_50","p1_95")
    return(res)
  }
  
  if(tm$uu<5 | tm$zz<5 ) { #| tm$uz<2
    nbin=30*nbin;eps=eps*1e-7; tolerance=10*tolerance; 
    if(verbose) print(c("bigger nbin",nbin,tolerance))}
  if(tm$zz<5){ p01 = 1.; lw01=0; dp01=1;} else {
    p01=seq(eps,1-eps,length.out = nbin); 
    lw01 = sapply(p01, ldp01, tm=tm); dp01=p01[2]-p01[1]}
  if(tm$uu<5){ p10 = 1.; lw10=0; dp10=1;} else {
    p10=seq(eps,1-eps,length.out = nbin); 
    lw10 = sapply(p10, ldp10, tm=tm); dp10=p10[2]-p10[1]}
  m = outer(p01,p10,function(x,y) x/(x+y))
  w = outer(lw01,lw10, function(x,y) {exp(x+y)*dp01*dp10})
  if(verbose) {
    print(sum(w))
    print(w)
  }
  if(!abs(sum(w)-1)<tolerance){ 
    print(sum(w))
    print(tm)}
  assertthat::assert_that(abs(sum(w)-1)<tolerance)
  mc = m[w>cutoff]
  wc = w[w>cutoff]
  ro = order(mc)
  ss = cumsum(wc[ro])
  m05 = mc[ro][which.min(abs(ss-.05))]
  m50 = mc[ro][which.min(abs(ss-.50))]
  m95 = mc[ro][which.min(abs(ss-.95))]
  res = c(m05,m50,m95)
  assertthat::assert_that(length(res)==3)
  names(res) = c("p1_05","p1_50","p1_95")
  #print(res)
  res
}

get_postsummary = function(tr) {
  post_summary(tmat(tr,stride=10), nbin=2000, tolerance=1e-2) %>% 
    as.list %>% as_tibble
}


#### ========= END Bayesian Markov Model estimation ==========



data4 = as_tibble(data3) %>% 
  mutate(postsumm = map(traj, get_postsummary))  %>%   
  unnest(postsumm) %>% select(-traj) %>% 
  filter(full_mean > 0.1)

data4 %>% ggplot(aes(x=full_mean, y=p1_50-full_mean)) + geom_point() +
  geom_errorbar(aes(ymin=p1_05-full_mean, 
                    ymax=p1_95-full_mean), 
                width=.01,
                position=position_dodge(.9))


### ==================  DO WATER MEDIATED INTERACTIONS ===================

## read the amber data, split by comma each bridge
raw = read_lines("hb.4k.formated.out") 
ee=str_split(raw,",")

## this converts each bridge string into a table with the involved residues and the frame number
## assume ab residues start at 858

f2m = function(string, frame){
  vals = sapply(str_split(str_squish(string)," "),as.numeric)
  vals.ab = vals[vals >= 858]
  vals.int = vals[vals < 858]
  bits = as_tibble(expand.grid(vals.int,vals.ab)) %>% 
    add_column(frame=frame)
  bits
}

## apply to all bridges in a frame
## apply to all frames in the trajectories
f2m_frame = function(frame_num) {do.call(rbind, lapply(ee[[frame_num]], f2m, frame=frame_num))}
all_data = do.call(rbind, lapply(1:4000, f2m_frame))

## generates the trajectory for a given interaction
get_traj = function(var1,var2){
  present = unlist(all_data %>% filter(Var1==var1, Var2 == var2) %>% select(frame))
  traj = 0*1:4000
  traj[present]=1
  traj
}

## generates the posterior summary from the trajectory; returns a tibble
get_postsummary = function(tr) {
  post_summary(tmat(tr,stride=10), nbin=2000, tolerance=1e-2) %>% 
    as.list %>% as_tibble
}

## just map get_traj and get get_postsummary; unnest and drop the trajectories
uu = all_data %>% group_by(Var1,Var2) %>% 
  summarise(count = n()) %>% filter(count>30) %>% 
  mutate(traj = map2(Var1, Var2, get_traj)) %>%
  mutate(postsumm = map(traj, get_postsummary)) %>% 
  unnest(postsumm) %>% select(-traj)

## recover the proper labels by joining with the name tables;
## calculate full_mean from counts;
## filter above 10%
uu2 = uu %>% mutate(int.residue = Var1, ab.residue = Var2) %>% 
  left_join(ab_names ) %>%
  left_join(int_names ) %>%
  mutate(full_mean = count/4000) %>%
  mutate(interaction = "Water_bridge") %>%
  filter(full_mean > .10)



#### PUT TOGETHER THE DATA FROM THE OTHER INTERACTIONS...
justthesevars = function(data) {
  data %>% select(interaction, 
                  contains("label"), contains("residue"), contains("chain"),
                  p1_05, p1_50, p1_95, full_mean)
}

#select(interaction, ab.label, int.label,  p1_05, p1_50, p1_95, full_mean)
data_complete = rbind(
  data4 %>%  justthesevars() %>% as_tibble(),
  uu2 %>% ungroup() %>% justthesevars())



data_complete %>% ggplot(aes(x=full_mean, y=p1_50-full_mean,
                         color=factor(interaction))) + 
  geom_point() +
  geom_errorbar(aes(ymin=p1_05-full_mean, 
                    ymax=p1_95-full_mean), 
                width=.01,
                position=position_dodge(.9))

### --------------------------- plot ---------------------------------------


data_complete %>% 
  mutate(label.ab = fct_reorder2(ab.label, ab.chain, ab.residue,.desc = FALSE),
         label.int = fct_reorder2(int.label, int.chain, int.residue, .desc = FALSE)
  ) %>% 
  ggplot(aes(x=label.ab, y=p1_50, group=interaction, fill=interaction)) + 
  geom_col(position="dodge",size=.3) +
  facet_wrap(.~label.int, scales = 'free_x', nrow=3) + 
  #facet_grid(.~label.int, scales = 'free_x', space='free_x') + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  #scale_fill_brewer(palette="RdYlBu", drop=TRUE) +
  scale_fill_manual(values=c('#fdae6b', #apolar
                             '#756bb1', #aro  weak purple
                             '#54278f', #elec dark purple 
                             '#08519c', #hb  b3
                             '#3182bd', #hb  b2
                             '#6baed6'  #wb  b1
  ), drop=TRUE) +
  geom_errorbar(aes(ymin=p1_05, #block_mean-0*block_sd, 
                    ymax=p1_95), #block_mean+block_sd), 
                width=.1,
                position=position_dodge(.9))

ggsave("facet_int.pdf", width = 9, height = 6)



data_complete %>% 
  mutate(label.ab = fct_reorder2(ab.label, ab.chain, ab.residue,.desc = FALSE),
         label.int = fct_reorder2(int.label, int.chain, int.residue, .desc = FALSE)
  ) %>% 
  ggplot(aes(x=label.int, y=p1_50, group=interaction, fill=interaction)) + 
  geom_col(position="dodge",size=.3) +
  facet_wrap(.~label.ab, scales = 'free_x', nrow=3) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  #scale_fill_brewer(palette="RdYlBu", drop=TRUE) +
  scale_fill_manual(values=c('#fdae6b', #apolar
                             '#756bb1', #aro  weak purple
                             '#54278f', #elec dark purple 
                             '#08519c', #hb  b3
                             '#3182bd', #hb  b2
                             '#6baed6'  #wb  b1
  ), drop=TRUE) +
  geom_errorbar(aes(ymin=p1_05, #block_mean-0*block_sd, 
                    ymax=p1_95), #block_mean+block_sd), 
                width=.1,
                position=position_dodge(.9))

ggsave("facet_ab.pdf", width = 9, height = 6)

################### generate pymol commands to color interface


######################
simplifyinteractions = function(data) {data %>% 
    mutate(simple.interaction = recode(interaction,
                                       Apolar='apolar', Aro_E2F='polar', Elec_ProP='polar',
                                       Hbond_ProA='polar', Hbond_ProD='polar', Water_bridge='polar'))
}



pymol_command = function(resID, chain, interac, strength){
  col = list(apolar='hotpink', polar='skyblue')[interac]
  if(interac=='polar') { 
    p1=colorRamp(head(brewer.pal(n = 9, "BuPu"),-1)) } else {
      p1=colorRamp(head(brewer.pal(n = 9, "Oranges"),-1)) }
  colname = paste0("col",resID,chain)
  s1=paste0("set_color ",colname,", [", paste(p1(strength),collapse = ','),"]; ")
  paste0(s1," color ",colname,", resi ",resID," and chain ",chain)
}

pymol_showsticks = function(resID, chain, interac, strength){
  paste0("show sticks, (resi ",resID," and chain ",chain, ") and not name c+o+n")
}


pp1= data_complete %>% 
  simplifyinteractions() %>% 
  group_by(int.label, simple.interaction) %>% 
  summarise(sint = sum(p1_50), mint=max(p1_50))%>%
  group_by(int.label) %>%
  summarise(
    dominant.interaction = simple.interaction[which.max(sint)], 
    mint = max(mint) ) %>% 
  left_join(int_names) %>% 
  mutate(newchain = recode(int.chain, aIIb='A', b3='B', H='C', L='D')) %>%
  mutate(chain_char = as.character(newchain), 
         interaction_char = as.character(dominant.interaction)) %>%
  mutate(command = pmap_chr(list(int.resID, chain_char, 
                                 interaction_char, mint), pymol_command),
         command2 = pmap_chr(list(int.resID, chain_char, 
                                  interaction_char, mint), pymol_showsticks)
  )


pp2= data_complete %>% 
  simplifyinteractions() %>% 
  group_by(ab.label, simple.interaction) %>% 
  summarise(sint = sum(p1_50), mint=max(p1_50))%>%
  group_by(ab.label) %>%
  summarise(
    dominant.interaction = simple.interaction[which.max(sint)], 
    mint = max(mint) ) %>% 
  left_join(ab_names) %>% 
  mutate(newchain = recode(ab.chain, aIIb='A', b3='B', H='C', L='D')) %>%
  mutate(chain_char = as.character(newchain), 
         interaction_char = as.character(dominant.interaction)) %>%
  mutate(command = pmap_chr(list(ab.resID, chain_char, 
                                 interaction_char, mint), pymol_command),
         command2 = pmap_chr(list(ab.resID, chain_char, 
                                  interaction_char, mint), pymol_showsticks))


## color residues
print(data.frame(pp1['command']), row.names = FALSE)
print(data.frame(pp2['command']), row.names = FALSE)

### show sticks
print(data.frame(pp1['command2']), row.names = FALSE)
print(data.frame(pp2['command2']), row.names = FALSE)



