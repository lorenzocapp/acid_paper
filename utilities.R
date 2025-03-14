##############################
####### UTILITIES ############
##############################


envelope = function(posterior, true)
{
  nn = length(posterior$grid)
  lo = posterior$q2.5
  hi = posterior$q97.5
  result = sum(true < hi & true > lo)
  return( result / nn)
}


bias = function(posterior, true)
{
  nn = length(posterior$grid)
  med = posterior$q50
  result = sum( (med - true)/true )
  return(result / nn)
}

dev = function(posterior, true)
{
  nn = length(posterior$grid)
  med = posterior$q50
  result = sum( abs(med - true)/true )
  return(result / nn)
}

relwid = function(posterior, true)
{

  nn = length(posterior$grid)
  lo = posterior$q2.5
  hi = posterior$q97.5
  result = sum( (hi - lo)/true )
  return(result / nn)
}



