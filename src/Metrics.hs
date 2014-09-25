-- | Calculate various metrics/statistics
module Metrics where

import AgrestiCoull
import MPileup (by_major_allele)
import Count
import Statistics.Distribution
import Statistics.Distribution.ChiSquared
import Data.List (foldl1', tails)

-- | Calculate vector angle between allele frequencies.  This is 
--   similar to `dist`, but from 1 (equal) to 0 (orthogonal)
angle :: Counts -> Counts -> Double
angle c1' c2' = let
  (c1,c2) = (toList c1', toList c2')
  vnorm = sqrt . sum . map ((**2))
  in sum $ zipWith (*) (map (/vnorm c1) c1) (map (/vnorm c2) c2)

-- Calculate pairwise nucleotide diversities
-- Not adjusting for coverage
ppi_params :: [Counts] -> [[Double]]
ppi_params (c:cs) = map (\x -> nd2 c x) (c:cs) : ppi_params cs
ppi_params [] = []

-- calculate diversity within and between sample pairs
fst_params :: [Counts] -> [[(Double,Double)]]
fst_params (x:xs) = go (x:xs)
  where go (y:ys) = map (heteroz $ y) ys : go ys
        go [] = []
fst_params [] = []

heteroz :: Counts -> Counts -> (Double,Double)
heteroz c1 c2 = let
  c1s = fromIntegral $ covC c1
  c2s = fromIntegral $ covC c2
  total = c1s + c2s
  hz x = 1 - fromIntegral (sq (getA x::Int) + sq (getC x) + sq (getG x) + sq (getT x))/fromIntegral (sq $ covC x) where sq z = z*z
  h_tot = hz (c1 `ptAdd` c2)
  h_subs = (hz c1*c1s + hz c2*c2s)/total
  in if c1s == 0 || c2s == 0 || h_tot == 0.0 then (0,0) 
     else (h_tot,h_subs)

heteroz_ :: [Double] -> [Double] -> (Double, Double)
heteroz_ c1 c2 = let
  c1s = sum c1
  c2s = sum c2
  total = c1s + c2s
  hz :: [Double] -> Double
  hz xs' = let s = sum xs' in 1 - sum (map ((**2) . (/s)) xs')
  h_tot, h_subs :: Double

  h_tot = hz $ zipWith (+) c1 c2
  h_subs = (hz c1*c1s + hz c2*c2s)/total
  in if c1s == 0 || c2s == 0 || h_tot == 0 then (0,0) 
     else (h_tot,h_subs)

-- | Calculate F_ST.  Note that this is weighted by the number of sequences (coverage)
--   this is not what we want for sequencing data!
f_st :: [Counts] -> Double
f_st xs = let
  h_subs, weights :: [Double]
  h_tot = nd (ptSum xs)
  h_subs = map nd xs
  weights = let t = fromIntegral $ covC (ptSum xs) in [fromIntegral (covC x)/t | x <- xs]
  in if h_tot == 0 then 0.0 
     else (h_tot - sum (zipWith (*) h_subs weights)) / h_tot

-- | Calculate F_ST - equivalent to the above
f_st_ :: [Counts] -> Double
f_st_ cs = let
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

-- | Calculate nucleotide diversity, the probability that sampling
-- twice will give you two different results.  Should we correct by
-- a factor of c/(c-1) here?  Note this gets weird with e.g. allele count of 
-- 1 and 1 (nd=1, rather than 0.5)
nd :: Counts -> Double
nd cs = let fs = pi_freqs cs
        in 1-sum (zipWith (*) fs fs)

-- | Nucleotide diversity between - the probability of selecting
--   one from each will differ.
nd2 :: Counts -> Counts -> Double
nd2 cs1 cs2 = let fs1 = pi_freqs cs1
                  fs2 = pi_freqs cs2
              in 1 - sum (zipWith (*) fs1 fs2)

-- | Calculate Pi, the probability of getting different nucleotides by
-- sampling from two different populations.
pi_k :: [Counts] -> Double
pi_k cs = sum [nd2 x y | (x:x2:xs) <- tails cs, y <- (x2:xs)] * 2/ lcs*(lcs-1)
  where lcs = fromIntegral (length cs)

-- pi_k_WRONG is the probability that sampling once from each
-- population will not be all the same.  One weakness is that if one
-- population has fifty-fifty allele frequencies, the result is always
-- exactly 0.5.  I.e. it can't identify divergent allele frequencies
-- in that case.  Like Fst, this also is indifferent to the actual
-- counts, so reliability depends on coverage.
pi_k_WRONG :: [Counts] -> Double
pi_k_WRONG cs = let fs = map pi_freqs cs
                    c = fromIntegral $ covC $ ptSum $ cs
  in if c>1 then c/(c-1)*(1 - (sum $ foldl1' (zipWith (*)) fs)) else 0

pi_freqs :: Counts -> [Double]
pi_freqs x = let s = fromIntegral $ covC x
             in [fromIntegral (getA_ x)/s,fromIntegral (getC_ x)/s,fromIntegral (getG_ x)/s,fromIntegral (getT_ x)/s]

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
--   allele frequency spectra.  Output * or +, depending on significance for each allele.
conf :: Counts -> Counts -> String
conf x y = let
  s1 = covC x  -- don't count structural variants
  s2 = covC y
  in [overlap (getA x,s1-getA x) (getA y,s2-getA y)
     ,overlap (getC x,s1-getC x) (getC y,s2-getC y)
     ,overlap (getG x,s1-getG x) (getG y,s2-getG y)
     ,overlap (getT x,s1-getT x) (getT y,s2-getT y)
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
--   the combined distribution, with error.
-- FIXME: use an error threshold min dist between major allele frequency conf intervals
conf_all :: [Counts] -> String
conf_all cs' = let
  -- attempt to smooth errors by subtracting one, seems to work:
  -- m x = max (x-1) 0
  -- rm_err x = C (0 `addA` (m $ getA x) `addC` (m $ getC x) `addG` (m $ getG x) `addT` (m $ getT x)) []
  in concat ["\t"++x | x <- map (conf (ptSum cs')) cs']

dsconf_pairs :: Double -> [Counts] -> String
dsconf_pairs e cs = go $ by_major_allele cs
  where go (x:xs) = " " ++ map (ds x) xs ++ go xs
        go [] = ""
        ds x y = if delta_sigma 2.326 x y > e then '*' else if delta_sigma 1.65 x y > e then '+' else '.'                                                                                                           
-- | Calculate distance (in absolute numbers) between confidence intervals 
--   with the given z-score.  This is overly conservative, and it is better to
--   use wald or wald_p below (or some more complex method, like Newcomb).
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
  (bs,bf) = (sum (map fst xs), sum (map snd xs))
  pairs = [((s,f),(bs-s,bf-f)) | (s,f) <- xs ]
  in map (uncurry (delta_sigma sig)) pairs

-- as above, only using wald_p instead of agresti-coull
dsw_all :: Double -> [Counts] -> [Double]
dsw_all sig counts = let
  xs = by_major_allele counts
  (bs,bf) = (sum (map fst xs), sum (map snd xs))
  pairs = [((s,f),(bs-s,bf-f)) | (s,f) <- xs ]
  in map (uncurry (wald_p sig)) pairs

-- | Wald intervals with pseudocounts
--   See Agresti and Caffo, 2000.
wald_p :: Double -> (Int, Int) -> (Int, Int) -> Double
wald_p z (s1,f1) (s2,f2) = wald z (s1+1,f1+1) (s2+1,f2+1)

-- | Calculate lower bound of Wald confidence interval for difference between frequencies
wald :: Double -> (Int, Int) -> (Int, Int) -> Double
wald z (s1,f1) (s2,f2) = let
  -- estimated success frequencies
  (//) x y = fromIntegral x / fromIntegral y
  n1 = s1+f1
  n2 = s2+f2
  p1 = s1//n1
  p2 = s2//n2
  in abs (p1-p2) - z*sqrt ((s1*f1)//(n1*n1*n1)+(s2*f2)//(n2*n2*n2))

-- | Calculate distance between approximate distributions
-- in terms of their standard deviation.  Perhaps use binomial distribution directly?
-- This is a z-score, i.e. score of 2 means that the 95% CIs barely overlap.
ci_dist :: Double -> (Int,Int) -> (Int,Int) -> Double
ci_dist z (s1,f1) (s2,f2) =
  let (i1,j1) = confidenceInterval z s1 f1
      (i2,j2) = confidenceInterval z s2 f2
      mu1 = i1+j1 -- all values are times two (so it cancels out)
      mu2 = i2+j2
      sd1 = j1-i1
      sd2 = j2-i2
  in if sd1+sd2 == 0 then 0 else abs (mu2-mu1)/(sd1+sd2)
