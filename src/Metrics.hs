-- | Calculate various metrics/statistics
module Metrics where

import AgrestiCoull
import MPileup (Counts(..), toList, sumList)
import Statistics.Distribution
import Statistics.Distribution.ChiSquared

-- | Calculate vector angle between allele frequencies.  This is 
--   similar to `dist`, but from 1 (equal) to 0 (orthogonal)
angle :: Counts -> Counts -> Double
angle c1' c2' = let
  (c1,c2) = (toList c1', toList c2')
  vnorm = sqrt . sum . map ((**2))
  in sum $ zipWith (*) (map (/vnorm c1) c1) (map (/vnorm c2) c2)

-- | Calculate F_ST
f_st :: [Counts] -> Double
f_st cs = let
  cs' = map toList cs
  hz :: [Double] -> Double
  hz xs' = let s = sum xs'
           in 1 - sum (map ((**2) . (/s)) xs')
  h_tot :: Double
  h_tot = hz $ sumList $ cs'
  h_subs, weights :: [Double]
  h_subs = map hz cs'
  weights = let total = sum $ concat cs'
            in [sum c/total | c <- cs']
  in if h_tot == 0 then 0.0 -- no heterozygosity in the population!
     else (h_tot - sum (zipWith (*) h_subs weights)) / h_tot

-- | Calculate Pi (my version), the expected number of differences
-- between two random samples from the populations.
pi_k :: [Counts] -> Double
pi_k cs' = let
  cs = map toList cs'
  no_diff = sum $ foldr1 (zipWith (*)) cs
  all_pairs = product $ map sum cs
  in (all_pairs - no_diff) / all_pairs

-- | Use AgrestiCoull to calculate significant differences between
--   allele frequency spectra
conf :: Counts -> Counts -> String
conf (C a1 c1 g1 t1 _v1) (C a2 c2 g2 t2 _v2) = let
  s1 = a1+c1+g1+t1  -- don't count structural variants
  s2 = a2+c2+g2+t2
  in [overlap (a1,s1-a1) (a2,s2-a2)
     ,overlap (c1,s1-c1) (c2,s2-c2)
     ,overlap (g1,s1-g1) (g2,s2-g2)
     ,overlap (t1,s1-t1) (t2,s2-t2)
     ]

-- | Helper function for conf
overlap :: (Int,Int) -> (Int,Int) -> Char
overlap (succ1,fail1) (succ2,fail2) =
  let (i1,j1) = confidenceInterval 1.65 succ1 fail1
      (i2,j2) = confidenceInterval 1.65 succ2 fail2
  in if i2>=j1 || i1>=j2 then
       let (k1,l1) = confidenceInterval 2.326 succ1 fail1
           (k2,l2) = confidenceInterval 2.326 succ2 fail2
       in if k2>=l1 || k1>=l2 then '*' else '+'
     else '.'

-- | Calculate distance between approximate distributions
-- in terms of their standard deviation.  Perhaps use binomial distribution directly?
-- This is a z-score, i.e. score of 2 means that the 95% CIs barely overlap.
ci_dist :: (Int,Int) -> (Int,Int) -> Double
ci_dist (s1,f1) (s2,f2) =
  let (i1,j1) = confidenceInterval 1.0 s1 f1
      (i2,j2) = confidenceInterval 1.0 s2 f2
      mu1 = i1+j1 -- all values are times two (so it cancels out)
      mu2 = i2+j2
      sd1 = j1-i1
      sd2 = j2-i2
  in if sd1+sd2 == 0 then 0 else abs (mu2-mu1)/(sd1+sd2)

-- | Calculate distance (in absolute numbers) between confidence intervals 
--   with the given z-score
delta_sigma :: Double -> (Int,Int) -> (Int,Int) -> Double
delta_sigma z (s1,f1) (s2,f2) =
  let (i1,j1) = confidenceInterval 1.0 s1 f1
      (i2,j2) = confidenceInterval 1.0 s2 f2
      mu1 = i1+j1 -- all values are times two (so it cancels out)
      mu2 = i2+j2
      sd1 = j1-i1
      sd2 = j2-i2
  in abs (mu2-mu1) - z*(sd1+sd2)

-- should probably include a warning if more than 20% of cells < 5 expected or some such
pearsons_chi² :: [[Int]] -> Double
pearsons_chi² t = let
  cols   = map sum
  rows x = if all null x then []
           else sum (map head x) : rows (map tail x)
  exps   = [[ fromIntegral (r*c) / fromIntegral (sum $ rows t) | r <- rows t] | c <- cols t ]
  chi    = sum [ (fromIntegral a-b)^(2::Int)/b | (a,b) <- zip (concat t) (concat exps) ]
  df     = (length t-1)*(length (head t) - 1)
  in if any (==0) (cols t) || any (==0) (rows t) then 1.0 else complCumulative (chiSquared df) chi

{-
-- | calcuate normalized vector distance between frequency counts
dist :: Counts -> Counts -> Double
dist c1 c2 = let
  vnorm = sqrt . sum . map ((**2))
  d1 = vnorm $ toList c1
  d2 = vnorm $ toList c2
  in vnorm $ [ x/d1-y/d2 | (x,y) <- zip (toList c1) (toList c2)]
-}

-- | Use AgrestiCoull to calculate significant difference from
--   a combined distribution, with error:
conf_all :: [Counts] -> String
conf_all cs' = let
  cs = map toList cs' :: [[Double]]
  [a,c,g,t] = map round  $ sumList cs -- add error!
  in concat ["\t"++x | x <- map (conf (C a c g t [])) cs']

{-
pseudo :: Int -> Int ->Counts -> Counts
pseudo tr err (C a c g t vs) =
  let s = tr + (a+c+g+t) `div` err
  in C (a+s) (c+s) (g+s) (t+s) vs
-}

