setwd("~/Desktop/integrin_antibody_figures/mutations")

files <- list.files(path = "results.1mutation", pattern='^[A-Z][1-9]*') 
read_data = function(filename,path='results.1mutation'){
  dt = read.table(paste0(path,'/',filename), header=FALSE)
  names(dt) = col_names = c("resn",'num','ddg','error')
  dt %>% add_column(name=filename)
}

mutation_data0 = do.call(rbind, lapply(files, read_data)) %>%
  extract(name, c("chain", "resi","resname"), '([:upper:])([:digit:]+)-([:upper:]+).+') %>%
  mutate( chain1 = case_when(chain=='A' ~ 'aIIb', chain=='B' ~ 'b3')) %>%
  mutate( res = aa321(resname)) %>% 
  mutate( mutres = aa321(resn)) %>%
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
    filter(absddg>2.) %>%
    mutate(mutant = reorder(mutant, ddg))

nthreshold = 4
mut1 = mutation_data %>% group_by(label) %>% top_n(nthreshold,ddg) %>% filter(ddg >0)
mut2 = mutation_data %>% group_by(label) %>% top_n(-nthreshold,ddg)%>% filter(ddg <0)

mut_plot = rbind(mut1,mut2)

mut_plot %>% 
  ggplot(aes(x=mutant, y=ddg, fill=ddg)) +
  geom_col(position = "dodge") +
  facet_wrap(.~label,scale='free_x', nrow=3) +
  scale_fill_distiller(palette = "PiYG") + #PuOr
  theme(axis.text.x=element_text(angle=90, hjust=1))  +
  geom_errorbar(aes(ymin=ddg-error, 
                    ymax=ddg+error), 
                width=.3,
                position=position_dodge(.9))
ggsave("mutations.pdf", width = 9, height = 6)

