setwd('~/Documents/code/skill_demonstartion/shallow_water/1d/Serial_Version/')
mydata <- read.table('data_0.txt', sep=" ")
plot(mydata[,1],mydata[,2],type="l",col='blue',lwd=3,xlab='location',ylab='height')
