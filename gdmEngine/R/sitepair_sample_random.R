#'@title Random site-pair sampler
#'
#'@description Sample site-pairs randomly. 
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param n.pairs.target (integer) The number of site-pairs to select.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns Dataframe, site-pairs table, being first 6 columns of the GDM input table, with dissimilarities not calculated.
#'
#'@examples output = sitepair_sample_random(My.site.env.data, n.pairs.target=10000, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom matrixStats rowMins rowMaxs
#'
#'@export
sitepair_sample_random=function(site.env.data,
                                n.pairs.target,
                                output.folder = NULL,       
                                output.name = "site_pairs_data",  
                                verbose=FALSE)
{
  # create a table to catch the row indices for the pairs selected for modelling
  train.pairs<-matrix(c(-3,-2,-1,0), nrow=2, ncol=2)
  colnames(train.pairs)<-c("temp.i", "temp.j")
  # now loop until we have enough pairs randomly sampled
  for(i in 1:1000)
  {
    # obtain row indices for randomly selected of pair-sites
    vals <- sample.int(nrow(site.env.data), (n.pairs.target*2), replace=TRUE)
    # format the random sample of site indices into pairs
    ij.pairs<-cbind(vals[c(1:n.pairs.target)], vals[c((n.pairs.target+1):(n.pairs.target*2))])         
    # rearrange indices for each site pair so that the smallest site index comes first
    temp.i<-rowMins(ij.pairs)
    temp.j<-rowMaxs(ij.pairs)
    ij.pairs<-cbind(temp.i,temp.j)
    # omit duplicate site pairs within this sample (note: we wouldn't do this step under a bootstrapping approach)
    ij.pairs<-unique(ij.pairs)
    # and omit site pairs that have already been selected
    ij.temp<-rbind(train.pairs, ij.pairs)
    ij.temp<-unique(ij.temp)
    selected.ij.pairs<-ij.temp[c((nrow(train.pairs)+1):nrow(ij.temp)),]
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
  
  return(Pairs.table)
  
}# end sitepair.sample.random
##---------------------------------------------------------------------------------------##
