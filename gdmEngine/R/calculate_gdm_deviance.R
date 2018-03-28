#'@title Calculate GDM deviance for observed & predicted dissimilarities
#'
#'@description Calculate GDM deviance for observed & predicted dissimilarities. Can be used for assessing cross-validation data. Translated from the c++ function 'CalcGDMDevianceDouble()' in the file 'NNLS_Double.cpp' from the GDM R package.
#'
#'@param pU (float) A vector of predicted dissimilarity values, of same length as pY.
#'@param pY (float) A vector of observeded dissimilarity values, of same length as pU.
#'
#'@return A single value (float) being the deviance.
#'
#'@examples output = calculate_gdm_deviance(My.pU, My.pY)
#'
#'@export
calculate_gdm_deviance=function(pU, 
                                pY)
{
dTotal = 0
for(i in 1:length(pU))
  {
  # calculate the first part (t1)
  if(pU[i] == 0) 
    {
    t1 = pY[i]
    } # end if 
  else
    {
    if(pY[i] == 0)
      {
      t1 = 0
      } # end if
    else
      {
      t1 = pY[i] * log( pY[i] / pU[i] )
      } # end else
    } # end else
  # calculate the second part (t2)
  if(pU[i] == 1)
    {
    t2 = 1 - pY[i]
    }# end if
  else 
    {
    if(pY[i] == 1)
      {
      t2 = 1 - pY[i]
      }# end if
    else
      {
      t2 = ( 1.0 - pY[i] ) * log( ( 1.0 - pY[i] ) / ( 1.0 - pU[i] ) )
      }# end else
    }# end else
  # accumulate the running sum
  dTotal <- dTotal +  ( t1 + t2 )
  }# end for i.row
  dTotal <- dTotal * 2
  return(dTotal)
} # end calculate_gdm_deviance
