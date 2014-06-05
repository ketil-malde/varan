{- 
  ESIV - Expected Site Information Value from SNPs

  Information value is the log odds for each allele, times the probablity of observing it
  not using any prior (or prior is 50/50)

  abs (avg(p1,p2)*log(p1/p2) - (1-avg(p1,p2))*log((1-p1)/(1-p2))
-}


module ESIV where

import Count
import AgrestiCoull

-- given a z-score `z` for confidence interval, and a
-- minimum error rate `epsilon`, calculate the ESIV conservatively
-- using the conf interval boundaries as frequencies.
esiv :: Double -> Double -> Counts -> Counts -> Double
esiv z epsilon c1 c2 = let
  t1 = covC c1
  t2 = covC c2
  esiv1 (a1,b1) (a2,b2) 
    | j1+epsilon<i2 = esiv_score (j1+epsilon/2) (i2-epsilon/2)
    | i1>j2+epsilon = esiv_score (j2+epsilon/2) (i1-epsilon/2)
    | otherwise     = 0
        where
          (i1,j1) = confidenceInterval z a1 b1
          (i2,j2) = confidenceInterval z a2 b2
  in sum [esiv1 (x,t1-x) (y,t2-y) | (x,y) <- zip (toList c1) (toList c2)]

esiv_score :: Double -> Double -> Double
esiv_score p1 p2 = abs ((p1+p2)/2*logBase 2 (p1/p2))

