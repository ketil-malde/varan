-- Calculates p-values using a permutation test (not currently used)

module RandomSelect where

import System.Random
import Control.Monad.State
import Variants
import AgrestiCoull (confidenceInterval)

pval_max_count :: Int
pval_max_count = 1000

pval_acc :: Double
pval_acc = 0.01

pval :: RandomGen g => g -> ([Counts] -> Double) -> [Counts] -> (Double,Double)
pval g f cs = let thresh = f cs
                  xs = take pval_max_count $ rselect g cs
                  go suc tot (y:ys) = let (a,b) = confidenceInterval 1.0 suc (tot-suc) 
                                      in if b-a <= pval_acc then fromIntegral suc/fromIntegral tot 
                                         else go (suc + if f y >= thresh then 1 else 0) (tot+1) ys
                  go suc tot [] = fromIntegral suc/fromIntegral tot
              in if thresh > 0 then (thresh, go 0 0 xs) else (0,1)

type AlleleSample = (Int,Int,Int,Int)

-- generate an infinite stream of sampled allele distributions
rselect :: RandomGen g => g -> [Counts] -> [[Counts]]
rselect g cs = map fromLists $ fst $ runState sel (g,sum_al,pool_sizes)
  where (_tot,sum_al,pool_sizes) = count cs
        
fromLists :: [AlleleSample] -> [Counts]
fromLists = map f
  where f (a,b,c,d) = C a b c d []
        
sel :: RandomGen g => State (g,AlleleSample,[Int]) [[AlleleSample]]
sel = do
  a <- select1
  as <- sel
  return (a:as)
        
select1 :: RandomGen g => State (g,AlleleSample,[Int]) [AlleleSample]
select1 = do
  -- select randomly for each p_sz
  (g,s_al,p_sz) <- get
  let (g',g'') = split g
  put (g'',s_al,p_sz) -- restore state
  return $ pickNs g' p_sz s_al

-- pick a random sampling of alleles given pool sizes
pickNs :: RandomGen g => g -> [Int] -> AlleleSample -> [AlleleSample]
pickNs _ [] _als = [] -- verify that als is empty
pickNs g (p:ps) als = let
  (g1,g2) = split g
  in pickN' g1 als p : pickNs g2 ps als

-- sample alleles given randomgen, count, and sum alleles (i.e. distribution)
pickN' :: RandomGen g => g -> AlleleSample -> Int -> AlleleSample
pickN' g als = go (0,0,0,0) g
  where
        go acc _ 0 = acc
        go acc g1 cnt = 
          let (i,g2) = randomR (1,sum' als) g1
          in go (add i als acc) g2 (cnt-1)
        sum' (a,b,c,d) = a+b+c+d
        
-- count total alleles, number of each variant, and number in each count
count :: [Counts] -> (Int,AlleleSample,[Int])
count cs = let xs = sumList . map toList $ cs
               ys = map (sum . toList) cs
               tuple [a,b,c,d] = (a,b,c,d)
           in (sum xs,tuple xs,ys)

add :: Int -> AlleleSample -> AlleleSample -> AlleleSample
add i (a,b,c,d) (x,y,z,w)
  | i <= a = (x+1,y,z,w)
  | i <= a+b = (x,y+1,z,w)
  | i <= a+b+c = (x,y,z+1,w)               
  | i <= a+b+c+d = (x,y,z,w+1)               
  | otherwise = error "too high i"

{-
-- remove allele number x from input
pick :: Int -> AlleleSample -> AlleleSample
pick x (a:as) | x <= a = 1 : map (const 0) as
              | otherwise = 0 : pick (x-a) as
pick _ [] = error "pick ran out of input, x too large"
-}