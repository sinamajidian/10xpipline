
library(ggplot2)
#library(ggpubr)

#setwd('/Volumes/Sina/University/WUR/Report-WUR/27/5m')


hap10="hap10"


# x_axis="Molecule length"
# list1=c(30,100,150)
# snp_rate="mollen" 



list1=c(2,3,5,7,10)
x_axis="Molecule number per bead"
snp_rate="molc" 





addr=paste('/Volumes/Sina/University/WUR/Report-WUR/30/',snp_rate, sep ='')
setwd(addr)



data_rr=data.frame(Methods='',coverage='',rr_sd='',rr='',stringsAsFactors=FALSE)
for (cov in list1){
  file_list=c(paste('pure_sdhap_rr',cov,'.txt', sep =''), 
             paste('wcc_rr',cov,'.txt', sep =''), 
              paste(hap10,'_rr',cov,'.txt', sep =''))
  
  for (i_file in 1:length(file_list)){
    val=integer(5) 
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")
      for (i in 1:5) { val[i]=mean(na.omit(as.numeric(data_read[i,]))) }
    }
    print(val)
    rr_sd=c(sd(val))
    rr_mean=c(mean(val))
    data_rr <- rbind(data_rr, data.frame(Methods=substring(file,1,4),coverage=cov,rr_sd=rr_sd,rr=rr_mean))
  }
}
data_rr
data_rr=data_rr[-1,]
data_rr$rr_sd= as.numeric(data_rr$rr_sd)
data_rr$rr= as.numeric(data_rr$rr)
data_rr$coverage= as.numeric(data_rr$coverage)
data_rr$coverage=factor(data_rr$coverage, levels = list1)
data_rr$Methods=factor(data_rr$Methods, levels =  c("pure","sdha","spn_","wcc_","hap1"))
levels(data_rr$Methods)[levels(data_rr$Methods)=="wcc_"] <- "Hap++"
levels(data_rr$Methods)[levels(data_rr$Methods)=="hap1"] <- "hap10"
levels(data_rr$Methods)[levels(data_rr$Methods)=="pure"] <- "SDhaP"

rownames(data_rr)=c()

data_rr_org=data_rr
#color1 = c("lightblue","darkgreen", "blue4", "maroon3") #"yellow",
color1 <- c("#999999", "#E69F00", "#56B4E9","#009E73","#F0E442","#0072B2")
            
          #, , "#D55E00", "#CC79A7")
#title="Simulated triploid from a 5mb-region solanum tuberosum",


data_rr_new=data_rr
data_rr1=read.table('/Volumes/Sina/University/WUR/Report-WUR/30/new_plot/tet_1m01a.txt',sep='\t',stringsAsFactors = F)
data_rr_new[13:15,]=data_rr[10:12,]
data_rr_new[13:15,2:4]=data.matrix(data_rr1[8:10,2:4])
data_rr=data_rr_new


rr_plot <-  ggplot(data_rr, aes(x=coverage, y=rr, fill=Methods))+
  geom_bar(width=0.8, position=position_dodge(width=0.88), stat="identity")+ 
  labs(x = x_axis, y ="Reconstruction rate")+
  geom_errorbar(aes(ymin=rr-rr_sd, ymax=rr+rr_sd),position=position_dodge(.9), width=.5) +
  scale_fill_manual(values = color1) +
  theme(legend.position="top",legend.title = element_blank())+
  ylim(0,1)+
  theme(legend.background = element_rect(size=0.6, linetype="solid",colour ="darkgrey")) #fill="lightblue"
rr_plot


#../new_plot/r
name_rr= paste('10rr',snp_rate,'.jpeg', sep ='')
jpeg(name_rr)
print(rr_plot)
dev.off()







### block length Bar plot
data_len=data.frame(Methods='',coverage='',l_sd='',BlockLength='',stringsAsFactors=FALSE)
for (cov in list1){  # length
  file_list=c(paste('pure_sdhap_length',cov,'.txt', sep =''), 
              paste('wcc_length',cov,'.txt', sep ='') ,
              paste(hap10,'_length',cov,'.txt', sep =''))

  
  for (i_file in 1:length(file_list)){
    length_block_mean=1+integer(5)
    sum_length=1+integer(5)
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")
      for (i in 1:5) {
        #line_number=na.omit(as.numeric(data_read[i,]))
        #length_i=integer(5)
       # for (j in 1:length(line_number)-1){
      #    length_i[j]=line_number[j+1]-line_number[j]
       # }
        length_i=na.omit(as.numeric(data_read[i,]))
        # if (length(length_i[length_i==0])){
        #   length_i[length_i==0]=NA
        #   length_i=na.omit(length_i)
        # }
      length_block_mean[i]=mean(length_i)
      sum_length[i]=sum(length_i)
      }
    }
    l_sd=c(sd(length_block_mean))
    l_mean=c(mean(length_block_mean))
    data_len <- rbind(data_len, data.frame(Methods=substring(file,1,4),coverage=cov,l_sd=l_sd,BlockLength=l_mean))
  }
}

data_len=data_len[-1,]
data_len

data_len$l_sd= as.numeric(data_len$l_sd)
data_len$BlockLength= as.numeric(data_len$BlockLength)
data_len$coverage= as.numeric(data_len$coverage)
data_len$coverage=factor(data_len$coverage, levels = list1)
data_len$Methods=factor(data_len$Methods, levels = c("pure","sdha","spn_","wcc_","hap1"))
levels(data_len$Methods)[levels(data_len$Methods)=="wcc_"] <- "Hap++"
levels(data_len$Methods)[levels(data_len$Methods)=="hap1"] <- "hap10"
levels(data_len$Methods)[levels(data_len$Methods)=="pure"] <- "SDhaP"




rownames(data_len)=c()



data_len_new=data_len
data10=read.table('/Volumes/Sina/University/WUR/Report-WUR/30/new_plot/tet_1m01a.txt',sep='\t',stringsAsFactors = F)
data_len_new[13:15,]=data_len[13:15,]
data_len_new[13:15,2:4]=data.matrix(data10[18:20,2:4])
data_len=data_len_new



len_plot <-  ggplot(data_len, aes(x=coverage, y=BlockLength, fill=Methods))+
  geom_bar(width=0.8, position=position_dodge(width=0.88), stat="identity")+
  labs(x = x_axis, y ="Mean haplotype block length (Log scale)")+ 
  geom_errorbar(aes(ymin=BlockLength-l_sd, ymax=BlockLength+l_sd),position=position_dodge(.9), width=.5)+
  scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values = color1) +
  theme(legend.position="top",legend.title = element_blank())+
  theme(legend.background = element_rect(size=0.6, linetype="solid",colour ="darkgrey")) #fill="lightblue"


len_plot

name_len= paste('10len',snp_rate,'.jpeg', sep ='')
jpeg(name_len)
print(len_plot)
dev.off()

















########## ver, please run previouse code for block length


#file_list=c('sdhap_ver5.txt','wcc_ver5.txt','hap10_2_alel1_ver5.txt')
#for (i_file in 1:length(file_list)){




data_ver=data.frame(Methods='',coverage='',ver_sd='', VER='',stringsAsFactors=FALSE)
for (cov in list1){
  file_list_len=c(paste('pure_sdhap_length',cov,'.txt', sep =''),
              paste('wcc_length',cov,'.txt', sep =''),
              paste(hap10,'_length',cov,'.txt', sep =''))
  file_list_ver=c(paste('pure_sdhap_ver',cov,'.txt', sep =''),
                  paste('wcc_ver',cov,'.txt', sep =''),
                 paste(hap10,'_ver',cov,'.txt', sep =''))
  

  for (i_file in 1:length(file_list_len) ){ # 1:length(file_list_len)
    file_len=file_list_len[i_file]
    file_ver=file_list_ver[i_file]
    if (file.exists(file_len)  && file.exists(file_ver) ) {
      data_read_len=read.csv(file_len, header = FALSE, sep = ",")
      data_read_ver=read.csv(file_ver, header = FALSE, sep = ",")
      ver=integer(5)
      for (i in 1:5) {
        length_i=na.omit(as.numeric(data_read_len[i,]))
        if (length(length_i[length_i==0])){
          length_i[length_i==0]=NA
          length_i=na.omit(length_i)
        }
       # print(length(length_i)-length(ve_i))
        
        length_i=as.numeric(length_i)
        ve_i=as.numeric(na.omit(as.numeric(data_read_ver[i,])))
        
        ve_i[ve_i==Inf]=length_i[ve_i==Inf]
        ver[i]=mean(as.numeric(mapply("/",ve_i,length_i,SIMPLIFY = FALSE)))
      } 
      #print(ver)
      ver_sd_=c(sd(ver))
      ver_mean_=c(mean(ver))
      data_ver <- rbind(data_ver, data.frame(Methods=substring(file_ver,1,4),coverage=cov,ver_sd=ver_sd_,VER=ver_mean_))
  }
  }
}

  

  

data_ver=data_ver[-1,]
data_ver

data_ver$ver_sd= as.numeric(data_ver$ver_sd)
data_ver$VER= as.numeric(data_ver$VER)
data_ver$coverage= as.numeric(data_ver$coverage)
data_ver$coverage=factor(data_ver$coverage, levels = list1)
data_ver$Methods=factor(data_ver$Methods, levels = c("pure","sdha","spn_","wcc_","hap1"))
levels(data_ver$Methods)[levels(data_ver$Methods)=="wcc_"] <- "Hap++"
levels(data_ver$Methods)[levels(data_ver$Methods)=="hap1"] <- "hap10"
levels(data_ver$Methods)[levels(data_ver$Methods)=="pure"] <- "SDhaP"

#data_ver$VER=-log(data_ver$VER)#data_ver$ver_sd=-log(data_ver$ver_sd)



rownames(data_ver)=c()
data_ver_org=data_ver





data_ver_new=data_ver
data10=read.table('/Volumes/Sina/University/WUR/Report-WUR/30/new_plot/tet_1m01a.txt',sep='\t',stringsAsFactors = F)
data_ver_new[13:15,]=data_ver[10:12,]
data_ver_new[13:15,2:4]=data.matrix(data10[28:30,2:4])
data_ver=data_ver_new

  
ver_plot <-    ggplot(data_ver, aes(x=coverage, y=VER, fill=Methods))+
  geom_bar(width=0.8, position=position_dodge(width=0.88), stat="identity")+
  labs( x = x_axis, y ="Vector Error Rate")+
  geom_errorbar(aes(ymin=VER-ver_sd, ymax=VER+ver_sd),position=position_dodge(.9), width=.5)+
  scale_fill_manual(values = color1)+
  theme(legend.position="top",legend.title = element_blank())+
  theme(legend.background = element_rect(size=0.6, linetype="solid",colour ="darkgrey")) 

ver_plot
# 

name_ver= paste('10ver',snp_rate,'.jpeg', sep ='')
jpeg(name_ver)
print(ver_plot)
dev.off()


name_out= paste('all',snp_rate,'.txt', sep ='')
write.table(data_rr, name_out, append = FALSE, sep = "\t",row.names = FALSE ) #
write.table(data_len, name_out, append = TRUE, sep = "\t",row.names = FALSE ) #
write.table(data_ver, name_out, append = TRUE, sep = "\t",row.names = FALSE ) #







