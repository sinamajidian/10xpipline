
library(ggplot2)


setwd('/Volumes/Sina/University/WUR/Report-WUR/27/5m')

#file_list=c('sdhap_rr2.txt','wcc_rr2.txt','spn_rr2.txt','hap10_rr2.txt')
#paste('wcc_rr',cov,'.txt', sep ='')




data_rr=data.frame(method='',coverage='',rr_sd='',rr='',stringsAsFactors=FALSE)
for (cov in c(2,5,10)){
  file_list=c(paste('sdhap_rr',cov,'.txt', sep =''),
              paste('wcc_rr',cov,'.txt', sep =''),
              paste('hap10_rr',cov,'.txt', sep =''))
  
  for (i_file in 1:length(file_list)){
    val=integer(5) 
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")
      for (i in 1:5) { val[i]=mean(na.omit(as.numeric(data_read[i,]))) }
    }
    rr_sd=c(sd(val))
    rr_mean=c(mean(val))
    data_rr <- rbind(data_rr, data.frame(method=substring(file,1,4),coverage=cov,rr_sd=rr_sd,rr=rr_mean))
  }
}
data_rr=data_rr[-1,]
data_rr


data_rr$rr_sd= as.numeric(data_rr$rr_sd)
data_rr$rr= as.numeric(data_rr$rr)
data_rr$coverage= as.numeric(data_rr$coverage)
data_rr$coverage=factor(data_rr$coverage, levels = c("2","5","10"))
data_rr$method=factor(data_rr$method, levels = c("sdha",'wcc_',"hap1"))
levels(data_rr$method)[levels(data_rr$method)=="sdha"] <- "pre1+sdhap"
levels(data_rr$method)[levels(data_rr$method)=="wcc_"] <- "pre1+pre2+sdhap"
levels(data_rr$method)[levels(data_rr$method)=="hap1"] <- "hap10"



rownames(data_rr)=c()
data_rr


rr_plot <-  ggplot(data_rr, aes(x=coverage, y=rr, fill=method))+
  geom_bar(position=position_dodge(), stat="identity")+
  labs(title="Simulated triploid from a 5mb-region solanum tuberosum", x = "Coverage", y ="Reconstruction rate")+
  geom_errorbar(aes(ymin=rr-rr_sd, ymax=rr+rr_sd),position=position_dodge(.9), width=.5)

rr_plot



jpeg("rr_001.jpeg")
print(rr_plot)  
dev.off() 



#matrix(0, 1, 3);


### block length

data_len=data.frame(method='',coverage='',l_sd='',BlockLength='',stringsAsFactors=FALSE)
for (cov in c(2,5,10)){
  file_list=c(paste('sdhap_line',cov,'.txt', sep =''),
              paste('wcc_line',cov,'.txt', sep =''),
              paste('hap10_line',cov,'.txt', sep =''))
  
  for (i_file in 1:length(file_list)){

    length_block_mean=1+integer(5) #matrix(NA,5,dim(data)[2]) #matrix(0, 1, 3);
    sum_length=1+integer(5)
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")

      for (i in 1:5) {
        line_number=na.omit(as.numeric(data_read[i,]))
        length_i=integer(5) #length(l)-1
        for (j in 1:length(line_number)-1){
          length_i[j]=line_number[j+1]-line_number[j]
        }
        if (length(length_i[length_i==0])){
          length_i[length_i==0]=NA
          length_i=na.omit(length_i)
        }
      
      #print(length_i)
      length_block_mean[i]=mean(length_i)
      sum_length[i]=sum(length_i)
      }
    }
    #print(mean(sum_length))
    
    l_sd=c(sd(length_block_mean))
    l_mean=c(mean(length_block_mean))
    data_len <- rbind(data_len, data.frame(method=substring(file,1,4),coverage=cov,l_sd=l_sd,BlockLength=l_mean))
  }
  
}
data_len=data_len[-1,]

data_len

data_len$l_sd= as.numeric(data_len$l_sd)
data_len$BlockLength= as.numeric(data_len$BlockLength)
data_len$coverage= as.numeric(data_len$coverage)
data_len$coverage=factor(data_len$coverage, levels = c("2","5","10"))
data_len$method=factor(data_len$method, levels = c("sdha",'wcc_',"hap1"))
levels(data_len$method)[levels(data_len$method)=="sdha"] <- "pre1+sdhap"
levels(data_len$method)[levels(data_len$method)=="wcc_"] <- "pre1+pre2+sdhap"
levels(data_len$method)[levels(data_len$method)=="hap1"] <- "hap10"



len_plot <-  ggplot(data_len, aes(x=coverage, y=BlockLength, fill=method))+
  geom_bar(position=position_dodge(), stat="identity")+
  labs(title="Simulated triploid from a 5mb-region solanum tuberosum", x = "Coverage", y ="Mean haplotype block length (Log scale)")+
  geom_errorbar(aes(ymin=BlockLength-l_sd, ymax=BlockLength+l_sd),position=position_dodge(.9), width=.5)+
  scale_y_continuous(trans = 'log10')

len_plot

jpeg("len.jpeg")
print(len_plot)  
dev.off() 

p <- ggplot(data_len, aes(x=coverage, y=BlockLength))+ #, fill=method
  geom_violin()

p












data_len=data.frame(method='',coverage='',BlockLength='',stringsAsFactors=FALSE)
for (cov in c(2,5,10)){
  file_list=c(paste('sdhap_line',cov,'.txt', sep =''),
              paste('wcc_line',cov,'.txt', sep =''),
              paste('hap10_line',cov,'.txt', sep =''))
  #file_list=c(paste('sdhap_line',cov,'.txt', sep =''))
  
  for (i_file in 1:length(file_list)){
    
    #length_block=matrix(1, 5, 3);
    #sum_length=1+integer(5)
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")
      
      for (i in 1:5) {
        line_number=na.omit(as.numeric(data_read[i,]))
        length_i=integer(5) #length(l)-1
        for (j in 1:length(line_number)-1){
          length_i[j]=line_number[j+1]-line_number[j]
        }
        if (length(length_i[length_i==0])){
          length_i[length_i==0]=NA
          length_i=na.omit(length_i)
        }
        #print(list(length_i))
        #length_block_mean[i]=mean(length_i)
        #sum_length[i]=sum(length_i)
      } # i 1:5 
      for (k in 1:length(length_i) ){
        data_len <- rbind(data_len, data.frame(method=substring(file,1,4),coverage=cov,BlockLength=length_i[k]))
      }
    } # if file exists
    #print(mean(sum_length))
    
    #l_sd=c(sd(length_block_mean))
    #l_mean=c(mean(length_block_mean))
    #data_len <- rbind(data_len, data.frame(method=substring(file,1,4),coverage=cov,l_sd=l_sd,BlockLength=l_mean))
  }
  
}
data_len=data_len[-1,]

data_len
#data_len=data_len[data_len$coverage!=2,]









data_len=data.frame(method='',coverage='',l_sd='',BlockLength='',stringsAsFactors=FALSE)

cov=5
#for (cov in c(2,5,10)){
 

# file_list=c(paste('sdhap_line',cov,'.txt', sep =''),
 #             paste('wcc_line',cov,'.txt', sep =''),
  #            paste('hap10_line',cov,'.txt', sep =''))
file_list=c(paste('sdhap_line',cov,'.txt', sep =''))

i_file=1
  for (i_file in 1:length(file_list)){
    
    length_block=matrix(0, 5, 3)
    file=file_list[i_file]
    if (file.exists(file)) {
      data_read=read.csv(file, header = FALSE, sep = ",")
      
      for (i in 1:5) {
        line_number=na.omit(as.numeric(data_read[i,]))
        length_i=integer(5) #length(l)-1
        for (j in 1:length(line_number)-1){
          length_i[j]=line_number[j+1]-line_number[j]
        }
        if (length(length_i[length_i==0])){
          length_i[length_i==0]=NA
          length_i=na.omit(length_i)
        }
        
        #print(length_i)
        length_block[i,i_file]=length_i
        sum_length[i]=sum(length_i)
      }
    }
    #print(mean(sum_length))
    
    l_sd=c(sd(length_block_mean))
    l_mean=c(mean(length_block_mean))
    data_len <- rbind(data_len, data.frame(method=substring(file,1,4),coverage=cov,l_sd=l_sd,BlockLength=l_mean))
  }
  
}





########## ver
data_rr=data.frame(method='',coverage='',sd='',mean='',stringsAsFactors=FALSE)
file_list=c('sdhap_ver5.txt','wcc_ver5.txt','hap10_alel1_ver5.txt')
for (i_file in 1:length(file_list)){
  file=file_list[i_file]
  data_read=read.csv(file, header = FALSE, sep = ",")
  val=integer(5) #matrix(0, 1, 3);
  # for (i in 2*(1:5)) { c5[i/2]=mean(na.omit(as.numeric(data[i,]))) }
  for (i in 1:5) { 
    val[i]=mean(na.omit(as.numeric(data_read[i,]))) 
    
  }
  rr_sd=c(sd(val))
  rr_mean=c(mean(val))
  data_rr <- rbind(data_rr, data.frame(method=file,coverage=c(5),sd=rr_sd,mean=rr_mean))
}


data_rr=data_rr[-1,]
data_rr$sd= as.numeric(data_rr$sd)
data_rr$mean= as.numeric(data_rr$mean)

rr_plot <- ggplot(data_rr, aes(x=method, y=mean,fill=method)) + #aes(x=dose, y=len, fill=supp)
  geom_bar(position=position_dodge(), stat="identity",width=.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)+ #,position=position_dodge(.9)
  labs(title="Triploid, coverage=5 perH", x = "Methods", y ="Reconstruction rate")

jpeg("rr.jpeg")
print(rr_plot)  
dev.off() 






