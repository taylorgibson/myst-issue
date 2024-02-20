# cleanNAs
# this function replaces all NAs
MyReplace = function(data) {data %>% t %>% na.locf(.,,T) %>% na.locf %>% t}
#shadenorm
shadenorm = function(below=NULL, above=NULL, pcts = c(0.025,0.975), mu=0, sig=1, numpts = 500, color = "gray", dens = 40,
                     lines=FALSE,between=NULL,outside=NULL){
  if(is.null(between)){
    bothnull = is.null(below) & is.null(above)
    if(bothnull==TRUE){
      below = ifelse(is.null(below), qnorm(pcts[1],mu,sig), below)
      above = ifelse(is.null(above), qnorm(pcts[2],mu,sig), above)
    }
  }

  if(is.null(outside)==FALSE){
    if(is.numeric(outside)==FALSE){if(outside==TRUE){outside=qnorm(pcts,mu,sig)}}
    below = min(outside)
    above = max(outside)
  }

  lowlim = mu - 4*sig
  uplim  = mu + 4*sig

  x.grid = seq(lowlim,uplim, length= numpts)
  dens.all = dnorm(x.grid,mean=mu, sd = sig)

  if(lines==FALSE){
    plot(x.grid, dens.all, type="l", xlab="X", ylab="Density")
  }

  if(lines==TRUE){
    lines(x.grid,dens.all)
  }

  if(is.null(below)==FALSE){
    x.below    = x.grid[x.grid<below]
    dens.below = dens.all[x.grid<below]
    polygon(c(x.below,rev(x.below)),c(rep(0,length(x.below)),rev(dens.below)),col=color,density=dens)
  }

  if(is.null(above)==FALSE){
    x.above    = x.grid[x.grid>above]
    dens.above = dens.all[x.grid>above]
    polygon(c(x.above,rev(x.above)),c(rep(0,length(x.above)),rev(dens.above)),col=color,density=dens)
  }



  if(is.null(between)==FALSE){
    if(is.numeric(between)==FALSE){if(between==TRUE){between=qnorm(pcts,mu,sig)}}
    from = min(between)
    to   = max(between)

    x.between    = x.grid[x.grid>from&x.grid<to]
    dens.between = dens.all[x.grid>from&x.grid<to]
    polygon(c(x.between,rev(x.between)),c(rep(0,length(x.between)),rev(dens.between)),col=color,density=dens)
  }
}
#transform data
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  rzT
}
########################################################
# functions to jazz up the pairs() all pairwise scatterplots
# some useful plotting functions
# see documentation for "pairs" function
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}
#
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2],0,1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


#########################################################
# correlation "scan" functions
# x is a matrix or dataframe with lots varaibles
# y is a vector
# cov is a covariate

# compute correlation of y with each element of x
mycorr <- function(x,y){
	apply(x,2,"cor",y,use="complete.obs")
}

# compute residuals to adjust x wrt covariate
myresids <- function(x, cov){
	residuals(lm(x~cov,na.action=na.exclude))
}

# compute correlation of y with each element of x after adjusting for covariate cov
mycorr.adj <- function(x,y,cov){
	x.adj <- apply(x,2,"myresids",cov)
	y.adj <- myresids(y,cov)
	apply(x.adj,2,"cor",y.adj,use="complete.obs")
}

# a FWER adjustment permutation test for correlation scans
mycorr.permute <- function(x,y,n){
	max.cor <- NULL
	for(i in 1:n){
			max.cor <- c(max.cor, max(abs(mycorr(x,sample(y))),na.rm=TRUE))
	}
	max.cor
}

# a FWER adjustment permutation test for adjusted correlation scans
mycorr.adj.permute <- function(x,y,cov,n){
	x.adj <- apply(x,2,"myresids",cov)
	y.adj <- myresids(y,cov)
	max.cor <- NULL
	for(i in 1:n){
			max.cor <- c(max.cor, max(abs(mycorr(x.adj,sample(y.adj))),na.rm=TRUE))
	}
	max.cor
}

R2Z <- function(r) log((1+r)/(1-r))

## random impute
rand.impute <- function(a) {
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- sample (a.obs, n.missing, replace=TRUE)
  return (imputed)
}

random.impute.data.frame <- function(dat, cols) {
  nms <- names(dat)
  for(col in cols) {
    name <- paste(nms[col],".imputed", sep = "")
    dat[name] <- rand.impute(dat[,col])
  }
  dat
}
#
simple.eda <-
  function(x) {
    ### create a simple function to explore data
    op <- par(no.readonly = TRUE); # save old parameters
    par(mfrow=c(1,3))
    hist(x);rug(x)
    boxplot(x);rug(x,side=2);title("boxplot")
    qqnorm(x);qqline(x)
    par(op)
  }