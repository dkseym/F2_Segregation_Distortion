#Get allele frequency estimates for each population
#2015-11-20 - Danelle Seymour
#Note - this code will need to be modified for specific uses (i.e. folder structure is specific to DKS)



# Beta-Binomial function
# In order to obtain the estimates for alpha and beta, the log-likelihood function needs to be defined.  Consult Bolker (In Press) for further information regarding the determination of the likelihood.

lbetabin = function (data, inits) {
	n <- data[,1] # This corresponds to the group_size
	y <-data[,2] # This corresponds to the data on incidence
	alpha <- inits[1] # This corresponds to the initial starting parameter for alpha
	beta <- inits[2] # This corresponds to the initial starting parameter for beta
	
	# Because optim minimizes a function, the negative log-likelihood is used. Also, the factorial is not necessary, as it does not depend on the parameters of interest
	sum(-lgamma(alpha+y)-lgamma(beta+n-y)+lgamma(alpha+beta+n)+lgamma(alpha)+lgamma(beta)-lgamma(alpha+beta))
}



#Create function for beta-binomial model of segregation distortion

setwd("./") # this directory structure will work if F2 allele frequency data is downloaded from the github respository (https://github.com/dkseym/F2_Segregation_Distortion)
files <- list.files("./allele_freq")

#GET CHROMOSOME LENGTHS
#>Chr1 CHROMOSOME dumped from ADB: Jun/20/09 14:53; last updated: 2009-02-02	30427671
#>Chr2 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02	19698289
#>Chr3 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02	23459830
#>Chr4 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02	18585056
#>Chr5 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02	26975502

chrlen <- c(30427671,19698289,23459830,18585056,26975502)
names(chrlen) <- c("Chr1","Chr2","Chr3","Chr4","Chr5")


#MODEL
Fit_BB<-function(input,window,step) {
	
	#Get data
	pos <- split(input$pos,input$chr)
	n <- split(input$n,input$chr)
	y <- split(input$y,input$chr)
	
	#Results storage
	err <- list()
	
	chr <- list()
	win <- list()
	
	alpha <- list()
	beta <- list()
	
	fitted.v <- list()
	upper.ci <- list()
	lower.ci <- list()
	
	pval <- list()
	
	#Loop through chromosomes
	for (c in names(pos)) {
		print(paste("Working on chromosome: ",c,sep=""))
	 	
	 	p <- pos[[c]]
	 	dd <- cbind(n[[c]],y[[c]])
	 	w <- seq(0,chrlen[c],step)
	 	
	 	e <- numeric(length(w))
	 	a <- numeric(length(w))
	 	b <- numeric(length(w))
	 	
	 	f <- numeric(length(w))
	 	u <- numeric(length(w))
	 	l <- numeric(length(w))
	 	pv <- numeric(length(w))
	 	
	 	#Loop through windows and fit model
	 	for (i in w) {

	 		r <- (i-window/2):(i+window/2-1)
	 		ddt <- matrix(dd[which(p%in%r),],ncol=2)
			inits <- c(10,10)
	
			#Run model
			tt <- try(optim(inits, lbetabin, method='L-BFGS-B',data=ddt,control=list(maxit=50),lower=1,upper=150))
			
			if (inherits(tt,"try-error")) {
				e[w==i] <- 1
				a[w==i] <- NA
				b[w==i] <- NA
				f[w==i] <- NA
				u[w==i] <- NA
				l[w==i] <- NA
				pv[w==i] <- NA
			}
	 		
	 		else if (tt$convergence==0) {
				e[w==i] <- tt$convergence
				a[w==i] <- tt$par[1]
				b[w==i] <- tt$par[2]
				f[w==i] <- tt$par[1]/(tt$par[1]+tt$par[2])
				u[w==i] <- qbeta(1-0.025,shape1=a[w==i],shape2=b[w==i])
				l[w==i] <- qbeta(0.025,shape1=a[w==i],shape2=b[w==i])
				pv[w==i] <- ks.test(0.5,function(x) pbeta(x,a[w==i],b[w==i]))$p.value
				
			}
			
			else {
				e[w==i] <- tt$convergence
				a[w==i] <- NA
				b[w==i] <- NA
				f[w==i] <- NA
				u[w==i] <- NA
				l[w==i] <- NA
				pv[w==i] <- NA
			}
		}
		#Combine chromosomes
		err[[c]] <- e
		
		chr[[c]] <- rep(c,length(w))
		win[[c]] <- w
		
		alpha[[c]] <- a
		beta[[c]] <- b
	
		fitted.v[[c]] <- f
		upper.ci[[c]] <- u
		lower.ci[[c]] <- l
	
		pval[[c]] <- pv
 
	 }
	#Output results
	results <- list()
	results[["chr"]] <- unlist(chr)
	results[["window"]] <- unlist(win,use.names=FALSE)
	results[["alpha"]] <- unlist(alpha,use.names=FALSE)
	results[["beta"]] <- unlist(beta,use.names=FALSE)
	results[["fitted.values"]] <- unlist(fitted.v,use.names=FALSE)
	results[["upper.ci"]] <- unlist(upper.ci,use.names=FALSE)
	results[["lower.ci"]] <- unlist(lower.ci,use.names=FALSE)
	results[["pvalues"]] <- unlist(pval,use.names=FALSE)
	results[["error"]] <- unlist(err,use.names=FALSE)
	return(results)
 }





##LOOP THROUGH ALL POPULATIONS
results <- list()
results.names <- c("chr","window","alpha","beta","fitted.values","upper.ci","lower.ci","pvalues","error")
results[results.names] <- list(NULL)
popid <- vector()

for (f in files) {
	
	#IMPORT DATA
	print(f)
	data <- try(read.table(paste("./allele_freq/",f,sep=""),header=FALSE,sep="\t"))
	
	if (inherits(data,"try-error")) {next}
	
	data[,10] <- round(data[,7]/data[,8])
	popid <- c(popid,strsplit(f,"_")[[1]][1])
	
	#FILTER DATA
	#REMOVE LINKED SNPS
	intersnp <- split(data[,3],data[,2])
	intersnp <- unlist(lapply(intersnp,diff))
	data <- data[-which(intersnp<100),]
	
	#OUTLIER REMOVAL
	quant <- quantile(data[,8])
	out.min <- quant[2]-1.5*(quant[4]-quant[2])
	out.max <- quant[4]+1.5*(quant[4]-quant[2])

	data <- data[which(data[,8]<out.max),]
	data <- data[which(data[,8]>out.min),]

	#RUN MODEL
	input <- data[,c(2,3,7,10)]
	colnames(input) <- c("chr","pos","y","n")
	
	out <- Fit_BB(input,5000000,500000)
	
	for (n in names(out)) {
		results[[n]] <- cbind(results[[n]],out[[n]])
	}
}


for (n in names(results)) {
	colnames(results[[n]]) <- popid
	write.table(results[[n]],paste("2015-11-20_model_bb_frequency_",n,".txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}
	


