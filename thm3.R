library("googlesheets4")
# googlesheets4::gs4_deauth()
# googlesheets4::gs4_auth()
# gs4_auth()

MINX=0
MAXX=20
SAMPLE_INTERVAL=1

genPlot<-function(r){
  rtext=paste("r",r,sep="")
  data<-read_sheet("https://docs.google.com/spreadsheets/d/1N6YMjxQiQEvTJs_nvwuPI4gxseqN1T2BwlPfQHgqdV4/edit?usp=sharing", sheet = rtext)
  data<-data[data$alpha > 0.5, ]
  # 
  # # Authenticate with Google Sheets
  # # This will open a browser window for authentication
  # # Follow the instructions to authenticate your R session with Google Sheets
  # # You only need to do this once
  # gs4_auth()
  
  # Example data for the first set of lines
  means1 <- data$mean_intersection[MINX:MAXX][seq(1, length(data$mean_intersection[MINX:MAXX]), by = SAMPLE_INTERVAL)]
  sds1 <- data$sd[MINX:MAXX][seq(1, length(data$sd[MINX:MAXX]), by = SAMPLE_INTERVAL)]
  xval <- data$alpha[MINX:MAXX][seq(1, length(data$alpha[MINX:MAXX]), by = SAMPLE_INTERVAL)]
  
  # Example data for the line plot
  line_data <- data$estimate[MINX:MAXX][seq(1, length(data$estimate[MINX:MAXX]), by = SAMPLE_INTERVAL)]
  
  
  
  titletext=paste(rtext)
  # Plotting the lines with error bars
  # plot(means1, ylim=c(min(c(means1 - 2*sds1, line_data)), max(c(means1 + 2*sds1, line_data))), 
  #      main=titletext,
  #      xlab="alpha", ylab="E|I| and L(1-q)", 
  #      type="l", col="blue")
  
  
  lines(xval, means1, ylim=c(min(c(means1 - 2*sds1, line_data)), max(c(means1 + 2*sds1, line_data))), type="l", col="blue")
  
  # Adding error bars for the first set
  for (i in 1:length(means1)) {
    arrows(xval[i], means1[i] - sds1[i], xval[i], means1[i] + sds1[i], angle=90, code=3, length=0.1, col="blue")
  }
  
  # Adding error bars for the second set
  # for (i in 1:length(means2)) {
  #   arrows(i, means2[i] - sds2[i], i, means2[i] + sds2[i], angle=90, code=3, length=0.1)
  # }
  
  # Adding the line plot alongside the error bars
  # lines(line_data, col="red")
  # Adding points to the line plot
  points(xval, line_data, col="red", pch=19)
  
  text(x=0, y=0.5, labels="L(1-q)", col="red")
  text(x=0, y=500, labels="E(|I|)", col="blue")
  
  # 
  # custom_labels <- data$k[0:10]
  # text(label_xs, par("usr")[3] - 0.5, labels = custom_labels, xpd = TRUE, srt = 90, adj = c(5, 0.5), col = "violet")
  
}

# Create an empty plot
plot(1, type = "n", xlim = c(MINX/2, MAXX/2), ylim = c(0, 10000), xlab = "alpha", ylab = "E(|I|) and L(1-q)", main = "r=0.01,0.05,0.2,0.5,L=10000", xaxt = "n", yaxt = "n")

axis(side = 1, at = seq(MINX/2, MAXX/2, by =SAMPLE_INTERVAL))


# xlab="alpha", ylab="E|I| and L(1-q)", 

r=0.01
genPlot(r)
# 
# # 
r=0.05
genPlot(r)


r=0.2
genPlot(r)

r=0.5
genPlot(r)