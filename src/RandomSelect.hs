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

type AlleleSample = [Int]

-- generate an infinite stream of sampled allele distributions
rselect :: RandomGen g => g -> [Counts] -> [[Counts]]
rselect g cs = map fromLists $ fst $ runState sel (g,sum_al,pool_sizes)
  where (_tot,sum_al,pool_sizes) = count cs
        
fromLists :: [AlleleSample] -> [Counts]
fromLists = map f
  where f [a,b,c,d] = C a b c d []
        f _ = error "Incorrect allelesample"
        
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
              