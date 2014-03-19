module MPileup (Counts(..), readPile1, toList, major_allele, by_major_allele, showC, showV, sumList) where

import Data.Char (toUpper)
import Data.List (foldl',intercalate,nub,elemIndex)

import Variants

-- convert counts to major/non-major allele counts
by_major_allele :: [Counts] -> [[Int]] -- always length 2
by_major_allele cs = let
  ls = map toList cs
  s  = sumList ls
  Just i = elemIndex (maximum s) s
  in map (\l -> [l!!i,sum l-l!!i]) ls

-- pick out major allele in first count, and output number of same/different
major_allele :: Counts -> Counts -> ((Int,Int),(Int,Int))
major_allele (C a b c d _) (C e f g h _) =
  let s1 = [a,b,c,d] 
      s2 = [e,f,g,h]
      m = maximum s1
      Just i = elemIndex m s1
  in ((m, sum s1-m),(s2!!i,sum [s2!!j | j <- [0..3], j /= i]))

-- count major allele in first sample
-- return flag whether informative, chrom, pos, ref, 
readPile1 :: String -> (Bool,String,String,Char,[Counts])
readPile1 = parse1 . words
  where
    parse1 (chr:pos:(ref:_):rest) = let trs = map snd $ triples ref rest
                                    in (check trs, chr,pos,ref,map count trs)
    parse1 xs = error ("parse1: insufficiently long line:"++show xs)
    
    check ts = let t = concat ts in null t || all (==head t) t

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


-- | Show SNP counts and coverage
showC :: Counts -> (String,Int)
showC (C as cs gs ts _) = (" "++(intercalate ":" $ map show [as,cs,gs,ts]),as+cs+gs+ts)

-- | Show structural variant count
showV :: [Counts] -> String
showV cs = let
  getv (C _ _ _ _ v) = v
  vs = nub $ concatMap getv cs
  countV :: Variant -> Counts -> Int
  countV v c = length . filter (==v) $ getv c
  in intercalate "\t" (show vs:[unwords $ map (\v -> show $ countV v c) vs | c <- cs])
