module MPileup (readPile, showPile, f_st) where

import Data.Char (toUpper)
import Data.List (foldl',intercalate,nub,elemIndex)

import AgrestiCoull
import Variants
import Text.Printf
import RandomSelect
import System.Random
import Control.Parallel

showPile :: (String,String,Char,[Counts]) -> IO String
showPile (_,_,_,[]) = error "Pileup with no data?"
showPile (chr,pos,ref,stats@(s1:ss)) = do
  g <- newStdGen
  let (f,pf) = pval g f_st (s1:ss)
      (p,pp) = pval g pi_k (s1:ss)
  pf `par` pp `pseq` return (
    chr++"\t"++pos++"\t"++[ref] ++concat ["\t"++showC s | s <- stats]
    ++"\t-"++concat ["\t"++conf s1 s | s <- ss] 
    --  ++ conf_all (s1:ss)
    ++concat [printf "\t%.3f" (angle s1 s) | s <- ss]
    ++print_pval (f,pf)
    ++print_pval (p,pp)
    ++concat [printf "\t%.2f" (uncurry ci_dist $ major_allele s1 s) | s <- ss]
    ++"\t"++showV stats)

-- pick out major allele in first count, and output number of same/different
major_allele :: Counts -> Counts -> ((Int,Int),(Int,Int))
major_allele (C a b c d _) (C e f g h _) =
  let s1 = [a,b,c,d] 
      s2 = [e,f,g,h]
      m = maximum s1
      Just i = elemIndex m s1
  in ((m, sum s1-m),(s2!!i,sum [s2!!j | j <- [0..3], j /= i]))

print_pval :: (Double, Double) -> String
print_pval (a,b) = printf "\t%.3f p=%.3f" a b

-- | calcuate normalized vector distance between frequency counts
dist :: Counts -> Counts -> Double
dist c1 c2 = let
  vnorm = sqrt . sum . map ((**2))
  d1 = vnorm $ toList c1
  d2 = vnorm $ toList c2
  in vnorm $ [ x/d1-y/d2 | (x,y) <- zip (toList c1) (toList c2)]

-- | similar to dist, but from 1 (equal) to 0 (orthogonal)
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

-- | Use AgrestiCoull to calculate significant difference from 
--   a combined distribution, with error:
conf_all :: [Counts] -> String
conf_all cs' = let
  cs = map toList cs'
  [a,c,g,t] = map round  $ sumList cs -- add error!
  p = id -- pseudo 1 30
  in concat ["\t"++x | x <- map (conf (p $ C a c g t [])) cs']

pseudo :: Int -> Int ->Counts -> Counts
pseudo tr err (C a c g t vs) = 
  let s = tr + (a+c+g+t) `div` err 
  in C (a+s) (c+s) (g+s) (t+s) vs

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

overlap :: (Int,Int) -> (Int,Int) -> Char
overlap (succ1,fail1) (succ2,fail2) = 
  let (i1,j1) = confidenceInterval 1.65 succ1 fail1
      (i2,j2) = confidenceInterval 1.65 succ2 fail2
  in if i2>=j1 || i1>=j2 then 
       let (k1,l1) = confidenceInterval 2.326 succ1 fail1
           (k2,l2) = confidenceInterval 2.326 succ2 fail2     
       in if k2>=l1 || k1>=l2 then '*' else '+'
     else '.'

-- calculate distance between approximate distributions
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

-- count major allele in first sample
-- return chrom, pos, ref, 
readPile :: String -> [(String,String,Char,[Counts])]
readPile = map (parse1 . words) . lines
  where
    parse1 (chr:pos:(ref:_):rest) = (chr,pos,ref,map (count . snd) $ triples ref rest)
    parse1 xs = error ("parse1: insufficiently long line:"++show xs)
    
    triples _ [] = []
    triples ref (cnt:bases:_quals:rest) = (cnt,parse ref $ map toUpper bases) : triples ref rest
    triples _ _ = error "triples: incorrect number of columns"
    
    count :: [Variant] -> Counts
    count = foldl' f (C 0 0 0 0 [])
      where f (C as cs gs ts vs) x = case x of
              Nuc 'A' -> (C (as+1) cs gs ts vs)
              Nuc 'C' -> (C as (cs+1) gs ts vs)
              Nuc 'G' -> (C as cs (gs+1) ts vs)
              Nuc 'T' -> (C as cs gs (ts+1) vs)
              Nuc 'N' -> (C as cs gs ts vs)
              Nuc _   -> error ("Not a nucleotide: "++show x)
              v -> (C as cs gs ts (v:vs))


-- | Show SNP counts
showC :: Counts -> String
showC (C as cs gs ts _) = " "++(intercalate " " $ map show [as,cs,gs,ts])++" "

-- | Show structural variant count
showV :: [Counts] -> String
showV cs = let
  getv (C _ _ _ _ v) = v
  vs = nub $ concatMap getv cs
  countV :: Variant -> Counts -> Int
  countV v c = length . filter (==v) $ getv c
  in intercalate "\t" (show vs:[unwords $ map (\v -> show $ countV v c) vs | c <- cs])
