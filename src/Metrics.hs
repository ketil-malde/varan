-- | Calculate various metrics/statistics
module Metrics where

import AgrestiCoull
import MPileup (Counts(..), toList, sumList, by_major_allele)
import Statistics.Distribution
import Statistics.Distribution.ChiSquared

-- | Calculate vector angle between allele frequencies.  This is 
--   similar to `dist`, but from 1 (equal) to 0 (orthogonal)
angle :: Counts -> Counts -> Double
angle c1' c2' = let
  (c1,c2) = (toList c1', toList c2')
  vnorm = sqrt . sum . map ((**2))
  in sum $ zipWith (*) (map (/vnorm c1) c1) (map (/vnorm c2) c2)

-- calculate total
fst_params :: [Counts] -> [[(Double,Double)]]
fst_params (x:xs) = go $ map toList (x:xs)
  where go (y:ys) = map (heteroz $ y) ys : go ys
        go [] = []
fst_params [] = []

heteroz :: [Double] -> [Double] -> (Double, Double)
heteroz c1 c2 = let
  total = sum (c1++c2)
  hz :: [Double] -> Double
  hz xs' = let s = sum xs'
           in 1 - sum (map ((**2) . (/s)) xs')
  h_tot :: Double
  h_tot = hz $ sumList $ [c1,c2]
  h_subs, weights :: [Double]
  h_subs = map hz [c1,c2]
  weights = [sum c/total | c <- [c1,c2]]
  in if total == 0 || sum c1 == 0 || sum c2 == 0 
     then (0,0) else (h_tot,sum (zipWith (*) h_subs weights))

-- | Calculate F_ST
f_st :: [Counts] -> Double
f_st cs = let
  cs' = map toList cs
  -- hm, er ikke dette bare nuc div?
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
-- between two random samples from the populations.  I.e. the
-- probability that sampling once from each population will not be all
-- the same.  One weakness is that if one population has fifty-fifty
-- allele frequencies, the result is always exactly 0.5.  I.e. it
-- can't identify divergent allele frequencies in that case.  Like Fst, this
-- also is indifferent to the actual counts, so reliability depends on coverage.
pi_k :: [Counts] -> Double 
pi_k cs = let fs = [ map (/sum x) x | x <- map toList cs]
  in 1 - (sum $ foldr1 (zipWith (*)) fs)

-- Or, equivalent
pi_k_alt :: [Counts] -> Double
pi_k_alt cs' = let
  cs = map toList cs'
  no_diff = sum $ foldr1 (zipWith (*)) cs
  all_pairs = product $ map sum cs
  in (all_pairs - no_diff) / all_pairs

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

-- | Use AgrestiCoull to calculate significant difference from
--   a combined distribution, with error.
conf_all :: [Counts] -> String
conf_all cs' = let
  cs = map toList cs' :: [[Double]]
  [a,c,g,t] = map round  $ sumList cs
  -- attempt to smooth errors by subtracting one, seems to work:
  m x = max (x-1) 0
  rm_err (C aa cc gg tt _) = (C (m aa) (m cc) (m gg) (m tt) [])
  in concat ["\t"++x | x <- map (conf (C a c g t []) . rm_err) cs']

-- | Calculate distance (in absolute numbers) between confidence intervals 
--   with the given z-score
delta_sigma :: Double -> (Int,Int) -> (Int,Int) -> Double
delta_sigma z (s1,f1) (s2,f2) =
  let (i1,j1) = confidenceInterval z s1 f1
      (i2,j2) = confidenceInterval z s2 f2
      mu1 = i1+j1 -- all values are times two (so it cancels out)
      mu2 = i2+j2
      sd1 = j1-i1
      sd2 = j2-i2
  in (abs (mu2-mu1) - (sd1+sd2))/2

-- use on output from by_major_allele
ds_all :: Double -> [Counts] -> [Double]
ds_all sig counts = let
  xs = by_major_allele counts
  (bs,bf) = (sum (map head xs), sum (map last xs))
  pairs = [((s,f),(bs-s,bf-f)) | [s,f] <- xs ]
  in map (uncurry (delta_sigma sig)) pairs

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
