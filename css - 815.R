rm(list=ls())
setwd('D:\\CSS\\CSS_数据')
data_list <- list.files(pattern = '.dta')
library(haven);library(tidyverse) 
# select variables
# generate a variable indicating year
css_list <- list()
lab_list_csv <- paste0(substr(data_list,1,7),'.csv')
lab_list <- list()
shiye_list <- list()
for (i in 1:length(data_list)) {
  css_list[[i]] <- read_dta(data_list[i]) %>% 
    rename_all(., .funs = tolower) 
  lab_list[[i]] <- css_list[[i]] %>% sjlabelled::get_label()
  shiye_list[[i]] <- lab_list[[i]][grep("失业",lab_list[[i]])]
  write.csv(lab_list[[i]],lab_list_csv[i],fileEncoding = 'gb2312')
}
# look for the variables labelled 'shiye'

# denoting the provinces with their names,not statistical codes.
shiye_list[[1]]

css2006 <- css_list[[1]] %>% select(qs2a,qc7a_06,qc7b_06,starts_with('qd13')) %>%
 mutate(province = labelled::to_factor(qs2a, levels = "labels") %>% 
          substr(1,2),
        wave = 2006,meet = 1*(qc7a_06 == 1),press = 1*(qc7a_06  %in% c(1,2)),
        soc_issue = ifelse(qd13_1==1|qd13_2==1|qd13_3==1,1,0)) %>% 
  select(province,wave,meet,press,soc_issue) 
table(css2006$meet == css2006$press) #meet and press is the same

css2006_sum <- css2006 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(press,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(press)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% ungroup()  %>% 
  column_to_rownames('province')

library(quest)
vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2006_sum), val = pat, pat = TRUE)})

css06 <- wide2long(css2006_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2006,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)


# Wave 2008
shiye_list[[2]]
css2008 <- css_list[[2]] %>% select(province,c9a7,c9b7,
                                    starts_with('e2')) %>% 
  mutate(meet = 1*(c9a7 == 1),
         press = 1*(c9b7 %in% c(1,2)),
         soc_issue = ifelse(e21==1|e22==1|e23==1,1,0) ) %>% 
  select(province,meet,press,soc_issue) 

css2008_sum <- css2008 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(press,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(press)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% ungroup()  %>% 
  column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
  str2str::pick(x = names(css2008_sum), val = pat, pat = TRUE)})


css08 <- wide2long(css2008_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2008,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)

# wave 2011
shiye_list[[3]]
css2011 <- css_list[[3]] %>% select(v41,f47,starts_with('g10')) %>% 
  mutate(province = substr(v41,1,2),meet = 1*(f47 == 1),
         soc_issue = ifelse(g101 == 1|g102 == 1|g103 == 1,1,0)) %>% 
  select(province,meet,soc_issue) 

css2011_sum <- css2011 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(soc_issue)),Sample = n()) %>% ungroup()  %>% 
  column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2011_sum), val = pat, pat = TRUE)})

css11 <- wide2long(css2011_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2011,Project = 'cgss', 
         Item = factor(obs,levels= 1:2,labels=c('meet','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)

# since 2013, the survey employ national code for province
# since 2013,there is a question: what prob to unemployment?
shiye_list[[4]]
apply(css_list[[4]] %>% select(b5,d48) ,2,table)

css2013 <- css_list[[4]] %>% select(prov_code,b5,d48,starts_with('f7')) %>% 
  mutate(province = prov_code,wave = 2013,meet = 1*(d48 == 1),
         unemp = 1*(b5 < 4),
         soc_issue = ifelse(f71 == 1|f72 == 1|f73 == 1,1,0)) %>% 
  select(province,wave,meet,unemp,soc_issue) 

library(rvest)
words <- read_html("https://www.converts.cn/citycode.html")
data <- html_table(html_nodes(words, "table"),
                   fill = TRUE,header=TRUE)[[1]] 

css2013$province <- factor(
  css2013$province,levels = data$行政区域代码,
  labels = substr(data$省级行政区,1,2) )
names(css2013)
css2013_sum <- css2013 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(unemp,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(unemp)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% ungroup()  %>% 
  column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2013_sum), 
                                      val = pat, pat = TRUE)})

css13 <- wide2long(css2013_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2013,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)

# wave 2015

shiye_list[[5]]  
css_list[[5]] %>% select(b5,d1_7,g5_1)
css2015 <- css_list[[5]] %>% select(n41,b5,d1_7,g5_1) %>% 
  mutate(province = substr(n41,1,2),wave = 2015,
         meet = 1*(d1_7 == 1),unemp = 1*(b5 < 4),
         soc_issue = ifelse(g5_1 == 1,1,0)) %>% 
  select(province,wave,meet,unemp,soc_issue) 

css2015_sum <- css2015 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(unemp,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(unemp)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% ungroup()  %>% 
  column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2015_sum), 
                                      val = pat, pat = TRUE)})

css15 <- wide2long(css2015_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2015,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)


# wave 2017

shiye_list[[6]]  
apply(css_list[[6]] %>% select(b5,d1_7,g5_1),2,table)  
css2017 <- css_list[[6]] %>% select(province,b5,d1_7,g5_1) %>% 
  mutate(wave = 2017,meet = 1*(d1_7 == 1),unemp = 1*(b5 < 4),
         soc_issue = ifelse(g5_1 == 1,1,0)) %>% 
  select(province,wave,meet,unemp,soc_issue) 
css2017[which(is.na(css_list[[6]]$province)),]
table(css_list[[6]]$province) %>% length() #two provinces is missing
prov_temp <- as.numeric(unlist(dimnames(table(css2017$province))))
code_temp <- data$省级行政区[which(data$行政区域代码 %in% prov_temp)]

css2017$province <- factor(
  css2017$province,levels = prov_temp,
  labels = substr(code_temp,1,2) ) 

mice::md.pattern(css2017)
css2017_sum <- css2017 %>% na.omit() %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(unemp,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(unemp)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% 
  ungroup()  %>% column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2017_sum), 
                                      val = pat, pat = TRUE)})

css17 <- wide2long(css2017_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2017,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)


# wave 2019
shiye_list[[7]]  
apply(css_list[[7]] %>% select(b5,d1_7,g5_1),2,table,useNA = "always") 
css_list[[7]]$prov_str %>% table()

css2019 <- css_list[[7]] %>% filter(psu != 1) %>% select(prov_str,b5,d1_7,g5_1) 
css2019[is.na(css2019)] = 0  # NA denotes zero

css2019 <- css2019 %>% mutate(province = substr(prov_str,1,2),wave = 2019,
                              meet = 1*(d1_7 == 1),unemp = 1*(b5 < 4),
                              soc_issue = ifelse(g5_1 == 1,1,0) ) %>% 
  select(province,wave,unemp,meet,soc_issue) 


css2019_sum <- css2019 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(unemp,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(unemp)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% 
  ungroup()  %>% column_to_rownames('province')



vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2019_sum), 
                                      val = pat, pat = TRUE)})

css19 <- wide2long(css2019_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2019,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)


# wave 2021
shiye_list[[8]] 

apply(css_list[[8]] %>% select(v3_1,b5g,d1a8,g4_1),2,table)
css2021 <- css_list[[8]] %>% select(v3_1,b5g,d1a8,g4_1) %>% 
  mutate(province = v3_1,wave = 2021,meet = 1*(d1a8 == 1),unemp = 1*(b5g < 4),
         soc_issue = ifelse(g4_1 == 1,1,0)) %>% 
  select(province,wave,meet,unemp,soc_issue) 

prov_temp <- as.numeric(unlist(dimnames(table(css2021$province))))
code_temp <- data$省级行政区[which(data$行政区域代码 %in% prov_temp)]

css2021$province <- factor(
  css2021$province,levels = prov_temp,
  labels = substr(code_temp,1,2) ) 



css2021_sum <- css2021 %>% group_by(province) %>% 
  summarise(Response.1 = mean(meet,na.rm = T),
            Response.2 = mean(unemp,na.rm = T),
            Response.3 = mean(soc_issue,na.rm = T),
            RespN.1 = sum(!is.na(meet)),
            RespN.2 = sum(!is.na(unemp)),
            RespN.3 = sum(!is.na(soc_issue)),Sample = n()) %>% 
  ungroup()  %>% column_to_rownames('province')

vrb_pat <- c("Response","RespN")
vrb_nm_list <- lapply(X = setNames(vrb_pat, nm = vrb_pat), 
                      FUN = function(pat) {
                        str2str::pick(x = names(css2021_sum), 
                                      val = pat, pat = TRUE)})

css21 <- wide2long(css2021_sum,vrb.nm = vrb_nm_list) %>% 
  rename(Country = Row.names) %>% 
  mutate(Year = 2021,Project = 'cgss', 
         Item = factor(obs,levels= 1:3,
                       labels=c('meet','press','soc_issue')),
         ItemCnt = Item) %>% 
  select(Country,Year,Project,Item,
         ItemCnt,Sample,Response,RespN)

df <- bind_rows(list(css06,css08,css11,css13,css15,css17,css19,css21)) %>% 
  arrange(Country,Year)





