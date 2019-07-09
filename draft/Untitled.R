
library(ggplot2)
library(dplyr)

setwd('/Volumes/Sina/University/WUR/Report-WUR/27/new')


########## rr
data_rr=data.frame(method='',coverage='',sd='',mean='',stringsAsFactors=FALSE)
file_list=c('sdhap_rr5.txt','wcc_rr5.txt','hap10_alel1_rr5.txt')
for (i_file in 1:length(file_list)){
  file=file_list[i_file]
  data_read=read.csv(file, header = FALSE, sep = ",")
  val=integer(5) #matrix(0, 1, 3);
  # for (i in 2*(1:5)) { c5[i/2]=mean(na.omit(as.numeric(data[i,]))) }
  for (i in 1:5) { val[i]=mean(na.omit(as.numeric(data_read[i,]))) }
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











### block length
block_length=matrix
data_len=data.frame(method='',coverage='',sd='',mean='',stringsAsFactors=FALSE)
file_list=c('sdhap_line5.txt','wcc_line5.txt','hap10_alel1_line5.txt')
for (i_file in 1:length(file_list)){
  file=file_list[i_file]
  data_read=read.csv(file, header = FALSE, sep = ",")
  length_block_mean=integer(5) #matrix(NA,5,dim(data)[2]) #matrix(0, 1, 3);
  for (i in 1:5) {
    line_number=na.omit(as.numeric(data_read[i,]))
    length_i=integer(5) #length(l)-1
    for (j in 1:length(line_number)-1){
      length_i[j]=line_number[j+1]-line_number[j]
    }
    length_block_mean[i]=mean(length_i)
  }
  len_sd=c(sd(length_block_mean))
  len_mean=c(mean(length_block_mean))
  data_len <- rbind(data_len, data.frame(method=file,coverage=c(5),sd=len_sd,mean=len_mean))
}


data_len=data_len[-1,]
data_len$sd= as.numeric(data_len$sd)
data_len$mean= as.numeric(data_len$mean)

len_plot<- ggplot(data_len, aes(x=method, y=mean,fill=method)) + #aes(x=dose, y=len, fill=supp)
  geom_bar(position=position_dodge(), stat="identity",width=.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)+ #,position=position_dodge(.9)
  labs(title="Triploid, coverage=5 perH", x = "Methods", y ="Haplotype block length")


jpeg("len.jpeg")
print(len_plot)  
dev.off() 





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






