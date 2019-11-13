# ----------------------
# ANALYSIS FOR BELL & PAEGLE
# Ethnic Markers and How to Find Them: An Ethnographic Investigation of Marker Presence, Recognition, and Social Influence
# University of Utah
# ----------------------
# packages you will need
require(jpeg)
require(graphics)
require(stats)
# graphics rendered on a mac
# so edit the graphics device
# call accordingly if on other OS


# data directory REPLACE TEXT WITH THE CORRECT DIRECTORY
datadir <- function(x="") paste0("Directory where the data files are","/",x)
setwd(datadir())
d <- read.csv( file=datadir("DesignsSurveyTonga2016short.csv"), header=TRUE, stringsAsFactors=FALSE, nrows=70)

# summary statistics of participants
sum( d$age >=1 & d$age<=2 )
numAgeClass <- function(x,l,u) {
	ans <- rep(NA,length(l))
	for( i in 1:length(l)) ans[i] <- sum( (x>=l[i] & x<=u[i]) )
	ans
	}
low <- c(18,31,51)
up <- c( 30, 50, 100)
numAgeClass( d$age, low,up)

# frequency of exposure by motif
freqY <- function(x) sum(x=="y")/(sum(x=="y")+sum(x=="n"))
freqY( d$des1seen )
freqY( d$des2seen )
seens <- paste0( "des",1:15,"seen")
named <- paste0( "des",1:15,"name")
freqsall <- sapply( 1:15, function(z) freqY( d[,seens[z]] ) )

namedmotif <- function( z ) sum( !d[,named[z] ]=="no" & !d[,named[z] ]=="" )/sum( !d[,named[z] ]=="" )

plotfreqs <- function( freqs, label=NULL ){
	or <- order(freqs, decreasing =TRUE)
	quartz(width=8.5, height=5)
	par( mar=c(0,10,1,1), las=1)	
	barplot( freqs[or], width=0.82, xlab="", ylim=c(-0.3,1), axes=FALSE, main=label )
	axis( side=2, at=c(0,1), las=1)
	mtext( "Fraction of sample\nrecognizing motif", side=2, las=1, line=9.5, at=0.5, adj=0)
	for ( i in 1:15 ){
		img<-readJPEG(datadir(paste0("m",or[i],".jpeg")))
		rasterImage(img,i-0.8,-0.1,i-0.2,-0.03)
		mtext( text=format( 100*namedmotif( or[i] ), digits=2), side=1, at=i-0.4, line=-3.5, adj=1, cex=0.7 )
		mtext( text="% named", side=1, line=-3.5, at=-1.5, cex=0.8)
		
	}
}

# FIGURE 2
plotfreqs( freqsall )

# ---------------------------------
# CLASSIFICATION
# ---------------------------------

# --------------------
# Classification tasks
# --------------------
# Load the motif classification data from Tonga
# loads object d.all.tonga 
load(datadir("tongaALLtriad2018.rdata"))
names(d.all.tonga)
names(d.all.tonga[[1]])
d.all.tonga[[1]]

# gender count
gc <- unlist(sapply( 1:length(d.all.tonga), function(z) d.all.tonga[[z]]$gender ))
sum( gc=="Tangata" )
sum( gc=="Fefine" )
length(d.all.tonga)

# age range
yob <- unlist(sapply( 1:length(d.all.tonga), function(z) d.all.tonga[[z]]$yob ))
# correct those that put age instead of year of birth, and large numbers
yob[82] <- 2018-26
yob[106] <- 2018-24
yob[114] <- 2018-20
yob[108] <- NA
yob[109] <- NA
ages <- 2018 - yob
summary(ages)

# Loads the classification data from Tonga
# loads object d.not
load(datadir("/utahtriad2018.rdata"))
names(d.not)
names(d.not[[1]])
# rename a nested data that matches
# the list names and structure of the Tonga data
for( i in 1:length(d.not) ){
	 d.not[[i]]$motif_triad <- d.not[[i]]$motif_triad_data[[1]]
	 d.not[[i]]$motif_tricomb <- d.not[[i]]$motif_tricomb[[1]]	 
	 } 
d.not[[1]]
# ---------------------------
# Calculate the round-by-round differences
# or average level of agreement

# function to count choice agreement
# between two triad tasks
# yields a number (Equation 1 in the manuscript)
triad.sim.count <- function( y1, t1, y2, t2){
	nObj <- 6
	tricomb <- t(combn(1:nObj,3))
	nr <- dim(tricomb)[1]
	y1tf <- ifelse( y1==1, TRUE, FALSE)
	y2tf <- ifelse( y2==1, TRUE, FALSE)	
	ans <- rep(NA,nr)
		for( i in 1:nr ){
			p1 <- tricomb[i,1]
			p2 <- tricomb[i,2]
			p3 <- tricomb[i,3]									
			add <- sapply( 1:nr, function(z) sum( sum(t1[z,] == p1), sum(t1[z,] == p2), sum(t1[z,] == p3 ) ) )
			rw <- which( add== 3 ) # which row
			choice1 <- t1[rw,y1tf[rw,2:4]]			
			choice2 <- t2[rw,y2tf[rw,2:4]]						
			ans[i] <- ifelse( choice1==choice2, 1, 0 )
		}	
	sum(ans)/20 # 20 total possible agreements
}



# function to organize the counting of triad task
# similarity counting between individuals
# yields an NxN matrix of counts
triad.sim <- function( d, dmn="motif" ){
	z <- paste0(dmn,"_triad")
	z1 <- paste0(dmn,"_tricomb")
	nObs <- length(d)	
	diff <- matrix( rep(NA,(nObs)^2), nrow=nObs )
	for( i in 1:nObs ){
		for( j in 1:nObs ) {
			y1 <- d[[i]][z][[1]]			
			triad1 <- d[[i]][z1][[1]]
			y2 <- d[[j]][z][[1]]		
			triad2 <- d[[j]][z1][[1]]
		if( !is.null(y1) & !is.null(triad1) & !is.null(y2) & !is.null(triad2) ) {			
			diff[i,j] <- triad.sim.count( y1 = y1, t1 = triad1, y2 = y2, t2 = triad2 )
			}else{
			diff[i,j] <- NA	
			}
		}		
	}
	diff 
}

# create similarity matrix bewteen US and Tonga participants
mtf.Tonga.US <- rep( list(0),sum( length(d.not), length(d.all.tonga ) ) )
for( i in 1:sum( length(d.not), length(d.all.tonga ) ) ){
 if( i<=length(d.not) ){ 
 	mtf.Tonga.US[[i]]$motif_triad <- d.not[[i]]$motif_triad
 	mtf.Tonga.US[[i]]$motif_tricomb <- d.not[[i]]$motif_tricomb
 	}else{
	mtf.Tonga.US[[i]]$motif_triad <- d.all.tonga[[i-length(d.not)]]$motif_triad
 	mtf.Tonga.US[[i]]$motif_tricomb <- d.all.tonga[[i-length(d.not)]]$motif_tricomb
	}
 }

length(mtf.Tonga.US)
mtf.Tonga.US[[1]]
mtf.tonga.us.mtrx <- triad.sim( d=mtf.Tonga.US, "motif" )

# Distance map
# Yes I know R has a built-in MDS function
# but I like to do it "by hand" just to check
DistanceMap = function( D, sub.names ){
	# Eigenvalues and eigenvectors
	decomp = eigen( D )
	decomp$values
	decomp$vectors

	# ----------------------------------------
	# MAP THE RELATIONSHIPS BETWEEN POPULATIONS
	# ----------------------------------------
	x.points = decomp$vectors[,1] * sqrt( decomp$values[1] ) # along first scaled eigenvector
	y.points = decomp$vectors[,2] * sqrt( decomp$values[2] ) # along second scaled eigenvector

	quartz(width=6, height=6)
	par( mar=c(3,3,3,3), bty="o" ) # no margins
	lim = max( x.points, y.points )
	xcenter <- median(x.points)
	xdist <-( max(x.points) - min(x.points) ) / 2
	ycenter <- median(y.points)	
	ydist <-( max(y.points) - min(y.points) ) / 2	
	plot( x.points, y.points, axes=F, type="n", ylab="", xlab="", xlim=c(xcenter-xdist,xcenter+xdist), ylim=c(ycenter-ydist,ycenter+ydist))
	axis( side=2, pos=xcenter )
	axis( side=1, pos=ycenter )
	points( x.points[sub.names=="T"], y.points[sub.names=="T"], lty=2, pch=1, xpd=TRUE)
	points( x.points[sub.names=="US"], y.points[sub.names=="US"], lty=2, pch=19, xpd=TRUE)	
	legend( "topright", box.lwd=0, pch=c(1,19), legend=c("Tonga","US"))
	}
	
mtxlabels <- c(rep("US",length(d.not)),rep("T",length(d.all.tonga) ) )	
DistanceMap( 1 - mtf.tonga.us.mtrx, sub.names <- mtxlabels  )


# Use other method for multidimensional scaling
loc <- cmdscale(1-mtf.tonga.us.mtrx)
x <- loc[, 1]
y <- -loc[, 2]

# FIGURE 3
plot(x, y, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE,
     main = "")
points(x, y, pch=ifelse(mtxlabels=="T", 1, 19), cex = 1)
legend( "topright", box.lwd=0.3, pch=c(1,19), legend=c("Tonga","US"))