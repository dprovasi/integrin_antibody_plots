setwd("~/Desktop/integrin_antibody_figures/mutations")
require(ggrepel)

aa321m = function (aa){
  aa1 <- c("-", ".", "X", bio3d::aa.table$aa1)
  aa3 <- c("---", "---", "UNK", bio3d::aa.table$aa3)
  aa1[aa3=='HIE']<-'He'
  aa1[aa3=='HID']<-'Hd'
  aa1[aa3=='HIP']<-'Hp'

  convert <- function(x) {
    if (is.na(x)) 
      return(NA)
    if (all(x != aa3)) {
      warning(paste("Unknown 3-letters code for aminoacid:", 
                    x))
      return("X")
    }
    else {
      return(aa1[which(x == aa3)])
    }
  }
  return(as.vector(unlist(sapply(aa, convert))))
}



files <- list.files(path = "results.1mutation", pattern='^[A-Z][1-9]*') 
read_data = function(filename,path='results.1mutation'){
  dt = read.table(paste0(path,'/',filename), header=FALSE)
  names(dt) = col_names = c("resn",'num','ddg','error')
  dt %>% add_column(name=filename)
}

mutation_data0 = do.call(rbind, lapply(files, read_data)) %>%
  extract(name, c("chain", "resi","resname"), '([:upper:])([:digit:]+)-([:upper:]+).+') %>%
  mutate( chain1 = case_when(chain=='A' ~ 'aIIb', chain=='B' ~ 'b3')) %>%
  mutate( res = aa321m(resname)) %>% 
  mutate( mutres = aa321m(resn)) %>%
  mutate( mutant = str_c(res, resi, mutres) ) %>%
  mutate(label= str_c(chain1, str_c(res, resi),sep = "/")) %>%
  mutate( absddg = abs(ddg)) %>%
  mutate( zscore = absddg/error) %>%
  mutate( pval = pt(zscore,df=20-1)) %>%
  drop_na() 

mutation_data0 %>% filter(ddg > 100) %>% select(label, mutant, ddg, error)


  mutation_data = mutation_data0 %>%
    filter(pval>0.7) %>%
    filter(absddg<20.) %>%
    filter(absddg>4) %>%
    mutate(mutant = reorder(mutant, ddg)) %>% as_tibble() %>%
    mutate(label = as_factor(label))
  
  
  ### as_factor and as.factor give different orders!!
  #%>%  mutate(label2 = fct_reorder2(factor(label), chain, resi,.desc = FALSE)) 
  #x=with(mutation_data, as_factor(label))
  #y=with(mutation_data,chain)
  #z=with(mutation_data,resi)
  #fct_reorder(x, z) 
  
  
nthreshold = 5
mut1 = mutation_data %>% group_by(label) %>% top_n(nthreshold,ddg) %>% filter(ddg >0)
mut2 = mutation_data %>% group_by(label) %>% top_n(-nthreshold,ddg)%>% filter(ddg <0)

mut_plot = rbind(mut1,mut2) 


#filter(label=="b3/M335" | label=="b3/D336") %>%
#  filter(label=="aIIb/E157") %>%
# filter(label=="b3/M335" | label=="b3/D336") %>%
  
mut_plot %>% ungroup() %>% 
  ggplot(aes(x=mutant, y=ddg, fill=ddg)) +
  geom_col(position = "dodge") +
  facet_wrap(.~label,scale='free_x',nrow=3) + ###  space='free' 
  scale_fill_distiller(palette = "PiYG") + #PuOr
  theme(axis.text.x=element_text(angle=90, hjust=1))  +
  geom_errorbar(aes(ymin=ddg-error, 
                    ymax=ddg+error), 
                width=.3,
                position=position_dodge(.9))
#ggsave("mutations-D336.pdf", width = 9, height = 6)
### facet_grid allows space=free


ggsave("mutations-marta.pdf", width = 9, height = 6)


###=========== needs pp1 from the FP

comparison = pp1 %>% select(-command, -command2) %>% 
  mutate(label=int.label) %>% 
  left_join(mut_plot ) 


comparison %>% mutate(
  chain = recode(int.chain, aIIb=21, b3=22),
  interaction_freq = mint) %>%
  ggplot(aes(x=ddg, y=interaction_freq, fill=ddg, label=mutant, shape=chain)) +
  geom_errorbarh(aes(xmin=ddg-error, 
                     xmax=ddg+error), 
                 height=.01, color='black',
                 position=position_dodge(.9)) +
  geom_point(aes(fill=ddg), colour='black', stroke=1, size=4) +
  scale_continuous_identity(aesthetics = 'shape', 
                            guide = 'legend', 
                            breaks=c(21,22),
                            labels=c("aIIb",'b3')) +
  scale_fill_distiller(palette = "PiYG") + #PuOr
  geom_text_repel(color='gray40',force=10) +
  ylim(0,1.2)
ggsave("comparison.pdf", width = 7, height = 6)


