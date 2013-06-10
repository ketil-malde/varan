module RandomSelect where

import System.Random
import Control.Monad.State
import Variants

type AlleleSample = [Int]

-- generate an infinite stream of sampled allele distributions
rselect :: RandomGen g => g -> [Counts] -> [[AlleleSample]]
rselect g cs = fst $ runState sel (g,tot,sum_al,pool_sizes)
  where (tot,sum_al,pool_sizes) = count cs
        
sel :: RandomGen g => State (g,Int,AlleleSample,[Int]) [[AlleleSample]]
sel = do
  a <- select1
  as <- sel
  return (a:as)
        
select1 :: RandomGen g => State (g,Int,AlleleSample,[Int]) [AlleleSample]
select1 = do
  -- select randomly for each p_sz
  (g,tot,s_al,p_sz) <- get
  let (g',g'') = split g
  put (g'',tot,s_al,p_sz) -- restore state
  return $ pickNs g' p_sz s_al

-- pick a random paritioning of alleles
pickNs :: RandomGen g => g -> [Int] -> AlleleSample -> [AlleleSample]
pickNs _ [] _als = [] -- verify that als is empty
pickNs g (p:ps) als = let
  (g1,g2) = split g
  (this,rest) = pickN' g1 p als
  in this : pickNs g2 ps rest
  
pickN' :: RandomGen g => g -> Int -> AlleleSample -> (AlleleSample,AlleleSample)
pickN' = go [0,0,0,0]
  where go acc _ 0 als = (acc,als)
        go acc g cnt als = 
          let (n,g') = randomR (1,sum als) g
              xs = pick n als
          in go (add xs acc) g' (cnt-1) (sub als xs)
        add = zipWith (+)
        sub = zipWith (-)

-- count total alleles, number of each variant, and number in each count
count :: [Counts] -> (Int,[Int],[Int])
count cs = let xs = sumList . map toList $ cs
               ys = map (sum . toList) cs
           in (sum xs,xs,ys)
        
-- remove allele number x from input
pick :: Int -> [Int] -> [Int]
pick x (a:as) | x <= a = 1 : map (const 0) as
              | otherwise = 0 : pick (x-a) as
pick _ [] = error "pick ran out of input, x too large"
              