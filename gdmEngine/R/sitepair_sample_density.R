#'@title Sample-density based site-pair sampler
#'
#'@description Sample site-pairs based on the number of nearby sites with composition data. 
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param n.pairs.target (integer) The number of site-pairs to select.
#'@param n.pairs.sample (integer) The number of pairs to assess simultaneously as part of the sampling process
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param a.used (float) The decay curve parameter (minimum y-value) for the number of times each site is used in selected pairs (default = 0.05)
#'@param b.used.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)
#'@param c.used (float) The decay curve parameter (slope of the curve) for the number of times each site is used in selected pairs. (default = 3) 
#'@param sigma.spair (float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = 0.5)
#'@param a.spair (float) The decay curve parameter (minimum y-value) for the number of sites within the neighbourhood radius (default = 0.05)
#'@param b.spair.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0
#'@param c.spair (float) The decay curve parameter (slope of the curve) for the number of sites within the neighbourhood radius (default = 3)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'site_pairs_data_dens'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe, site-pairs table, being first 6 columns of the GDM input table, with dissimilarities not calculated.
#'
#'@examples output = sitepair_sample_density(My.site.env.data, n.pairs.target=10000, domain.mask=My.mask)
#'
#'@importFrom matrixStats rowMins rowMaxs
#'@importFrom spatstat owin ppp density.ppp
#'
#'@export
sitepair_sample_density=function(site.env.data,
                                    n.pairs.target,
                                    n.pairs.sample=NULL, 
                                    domain.mask,
                                    a.used=0.05, 
                                    b.used.factor=2, 
                                    c.used=3, 
                                    sigma.spair=0.5,
                                    a.spair=0.05, 
                                    b.spair.factor=1.0, 
                                    c.spair=1, 
                                    output.folder = NULL,       
                                    output.name = "site_pairs_data_dens",  
                                    verbose=FALSE)
{
  ## Establish the working parameters ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  n.pairs.total <- ((nrow(site.env.data)^2)-nrow(site.env.data))/2
  if(is.null(n.pairs.sample))
    {
    n.pairs.sample<-floor(n.pairs.target/10)
    }# end if is.null(n.pairs.sample)
  b.used<-(n.pairs.target/nrow(site.env.data))*b.used.factor
  # Calculate the density of sites around each point, to use as a weighting
  sites.window <- owin(xrange=c(domain.mask@extent@xmin,domain.mask@extent@xmax),yrange=c(domain.mask@extent@ymin,domain.mask@extent@ymax))
  sites.points <- ppp(x=site.env.data$decimalLongitude, y=site.env.data$decimalLatitude, window=sites.window)
  sites.density <- density.ppp(sites.points, sigma=sigma.spair, at="points", leaveoneout=FALSE, positive=TRUE, verbose=TRUE)
#plot(x=site.env.data$decimalLongitude, y=site.env.data$decimalLatitude, cex=10/sites.density)
  b.spair<-mean(sites.density)*b.spair.factor
  SiteDensity.Wt<-decay.curve(sites.density, a.spair, b.spair, c.spair)
  # Create a table to catch the row indices for the pairs selected for modelling
  train.pairs<-matrix(c(-3,-2,-1,0), nrow=2, ncol=2)
  colnames(train.pairs)<-c("temp.i", "temp.j")
  # Set up the site weighting table (site.ID, Dist2NSW.Wt, ntimes.used, PairUse.Wt)
  train.plot.weight.table <- data.frame("site.ID" = site.env.data$xy, 
                                        "ntimes.used" = 0, 
                                        "PairUse.Wt" = decay.curve(0, a.used, b.used, c.used),
                                        "SiteDensity.Wt" = SiteDensity.Wt)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  ## Implement the sampling code, looping over sets of candidate pairs, sampling those,  and continuing until the
  ## specified number of pairs has been sampled (or it just takes too long - >1000 candidate samples). ~~~~~~~~##
  # now loop until we have enough pairs randomly sampled
  for(i in 1:1000)
    {
    # obtain row indices for randomly selected of pair-sites
    vals <- sample.int(nrow(site.env.data), (n.pairs.sample*2), replace=TRUE)
    # format the random sample of site indices into pairs
    ij.pairs<-cbind(vals[c(1:n.pairs.sample)], vals[c((n.pairs.sample+1):(n.pairs.sample*2))])         
    # rearrange indices for each site pair so that the smalles comes first
    temp.i<-rowMins(ij.pairs)
    temp.j<-rowMaxs(ij.pairs)
    ij.pairs<-cbind(temp.i,temp.j)
    # omit duplicate site pairs within this sample (note: we wouldn't do this step under a bootstrapping approach)
    ij.pairs<-unique(ij.pairs)
    # and omit site pairs that have already been selected
    ij.temp<-rbind(train.pairs, ij.pairs)
    ij.temp<-unique(ij.temp)
    ij.pairs<-ij.temp[c((nrow(train.pairs)+1):nrow(ij.temp)),]
    # note how many unique sample pairs we've selected
    n.pairs.selected<-nrow(ij.pairs)
    # And finally, calculate the total weight for each pair, combining the distance to NSW weights, 
    # the ntimes used weights, and the distance between sites in the pair weight
    Site.Pair.Prob <- train.plot.weight.table$PairUse.Wt[ij.pairs[,1]] * train.plot.weight.table$PairUse.Wt[ij.pairs[,2]] * train.plot.weight.table$SiteDensity.Wt[ij.pairs[,1]] * train.plot.weight.table$SiteDensity.Wt[ij.pairs[,2]]
    # Use these probabilities to randomly select site-pairs
    Pair.Sample <- rbinom(n=length(Site.Pair.Prob), size=1, prob=Site.Pair.Prob)
    selected.ij.pairs <- ij.pairs[Pair.Sample>0,]
    # update the weight's table ntimes selected for these sites
    i.freq <- table(selected.ij.pairs) 
    i.freq <- as.data.frame(i.freq)
    i.freq$selected.ij.pairs <- as.numeric(as.character(i.freq$selected.ij.pairs)) # 
    train.plot.weight.table$ntimes.used[i.freq$selected.ij.pairs] <- train.plot.weight.table$ntimes.used[i.freq$selected.ij.pairs] + i.freq$Freq
    train.plot.weight.table$PairUse.Wt <- decay.curve(train.plot.weight.table$ntimes.used, a.used, b.used, c.used)
    # and add these pairs to the main list of pairs selected for modelling
    train.pairs<-rbind(train.pairs, selected.ij.pairs)
    # check if we have enough pairs now (if so, randomly remove the necessary amount, so we hit our target,
    # then break out of the loop)
    if(nrow(train.pairs) >= (n.pairs.target + 2))
      {
      # remove the initial rows from train.pairs
      train.pairs<-train.pairs[-c(1,2),]
      # check how many excess pairs we have
      n.excess <- nrow(train.pairs) - n.pairs.target 
      # Randomly select pairs to drop & remove them from the selected list
      if(n.excess > 0)
      {
        drop.indices <- sample(seq_len(n.excess), size = n.excess, replace = FALSE)
        train.pairs<-train.pairs[-drop.indices,]
      }# end if n.excess > 0
      # And break out of the loop
      break()
      }# end if nrow(train.pairs) >= (n.pairs.target + 2)
    } # end for i
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  ## Processing the selected pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Provide a warning if not enough pairs were selected
  if(nrow(train.pairs) < n.pairs.target)
  {
    # remove the initial rows from train.pairs
    train.pairs<-train.pairs[-c(1,2),]
    warning("less pairs selected than desired", call. = FALSE)
  }
  # Prepare the start of a GDM input table for the pairs selected
  Pairs.table <- data.frame(distance	= 0,
                            weights = 1,
                            s1.xCoord = site.env.data$decimalLongitude[train.pairs[,1]],
                            s1.yCoord = site.env.data$decimalLatitude[train.pairs[,1]],
                            s2.xCoord = site.env.data$decimalLongitude[train.pairs[,2]],
                            s2.yCoord = site.env.data$decimalLatitude[train.pairs[,2]]) 
  # return the selected pairs
  return(Pairs.table)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
}# end sitepair_sample_geographic
##-------------------------------------------------------------------------------------------------------------##
