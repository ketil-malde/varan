module MPileup where

import Data.Char (toUpper)
import Data.List (foldl',intercalate,nub)

import AgrestiCoull
import Variants
import Text.Printf

showPile :: (String,String,Char,[Counts]) -> String
showPile (chr,pos,ref,stats@(s1:ss)) = 
  chr++"\t"++pos++"\t"++[ref] ++concat ["\t"++showC s | s <- stats]
  ++"\t-"++concat ["\t"++conf s1 s | s <- ss] 
  ++concat [printf "\t%.3f" (dist s1 s) | s <- ss]
  ++printf "\t%.3f" (f_st (s1:ss))
  ++"\t"++showV stats

-- | calcuate normalized vector distance between frequency counts
dist :: Counts -> Counts -> Double
dist (C a1 c1 g1 t1 _v1) (C a2 c2 g2 t2 _v2) = let
  vnorm = sqrt . sum . map ((**2))
  d1 = vnorm $ map fromIntegral [a1,c1,g1,t1]
  d2 = vnorm $ map fromIntegral [a2,c2,g2,t2]
  in vnorm $ [ fromIntegral x/d1-fromIntegral y/d2 | (x,y) <- zip [a1,c1,g1,t1] [a2,c2,g2,t2]]

f_st :: [Counts] -> Double
f_st cs = let
  toList (C x y z v _) = map fromIntegral [x,y,z,v]
  cs' = map toList cs
  hz :: [Double] -> Double
  hz xs' = let s = sum xs'
           in 1 - sum (map ((**2) . (/s)) xs')
  sumList = foldr (zipWith (+)) [0,0,0,0]            
  h_tot :: Double
  h_tot = hz $ sumList $ cs'
  h_subs, weights :: [Double]
  h_subs = map hz cs'
  weights = let total = sum $ concat cs'
            in [sum c/total | c <- cs']
  in (h_tot - sum (zipWith (*) h_subs weights)) / h_tot

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


data Counts = C !Int !Int !Int !Int [Variant]

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
