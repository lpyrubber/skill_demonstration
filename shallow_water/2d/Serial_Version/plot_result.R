current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
mydata <- read.table('data_0.txt', sep=' ')
x <- mydata[,1]
y <- mydata[,2]
z <- mydata[,3]
n <- sqrt(length(x))
xm <- min(x)
xM <- max(x)
x <- seq(xm, xM, length = n)
n <- sqrt(length(y))
xm <- min(y)
xM <- max(y)
y <- seq(xm, xM, length = n)
library(plotly)
z <- matrix(z,nrow = n, ncol = n)
fig <- plot_ly(x=x,y=y,z=z)
fig <- fig %>% add_surface
fig
