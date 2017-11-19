######## 
# GBIF #
########


install.packages('rgbif')
library(rgbif)
?occ_count

## Perhaps there is a 200 K limit on downloads? 
## Or some series handbreak for requests greater than this?
## See https://dagendresen.wordpress.com/2015/09/11/downloading-occurrence-data-using-the-gbif-rest-api/

## Try with class:Equisetopsida

name_lookup(query='Equisetopsida', rank = 'CLASS', limit=1)
# nubKey should relate to the GBIF taxonomic backbone. 
occ_count(taxonKey = 246) # ... < 1M records...

name_lookup(query='mammalia', limit=1)
occ_count(taxonKey = 359) # ~ 12 M records...

name_lookup(query='Plantae', rank = 'KINGDOM', limit=1)
# nubKey should relate to the GBIF taxonomic backbone. 
occ_count(taxonKey = 6) # ... lots: 205278309