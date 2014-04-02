{-# Language BangPatterns #-}
module MPileup (Counts(..), readPile1, toList, major_allele, by_major_allele, showC, showV, sumList, MPileRecord) where

import Data.Char (toUpper)
import Data.List (intercalate,nub,elemIndex)
import qualified Data.ByteString.Lazy.Char8 as B

import Variants hiding (parse)

type MPileRecord = (Bool,B.ByteString,B.ByteString,Char,[Counts])

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
readPile1 :: B.ByteString -> MPileRecord 
readPile1 = parse1 . B.words
  where
    parse1 (chr:pos:r:rest) = let trs = triples ref rest
                                  ref = B.head r
                              in (ignore trs, chr, pos, ref, trs)
    parse1 xs = error ("parse1: insufficiently long line:"++show xs)
    
    ignore cs = (<=1) $ length $ filter (/=(0::Int)) $ sumList $ map toList cs

    triples _ [] = []
    triples ref (_cnt:bases:_quals:rest) = parse ref (C 0 0 0 0 []) (B.map toUpper bases) : triples ref rest
    triples _ _ = error "triples: incorrect number of columns"
    
    -- this could probalby be faster if counting was incorporated directly, 
    -- avoiding the intermediate [Variant] data structure
    parse :: Char -> Counts ->  B.ByteString -> Counts
    parse ref !cts bs = case B.uncons bs of
      Nothing -> cts
      Just (c_,str_) -> p c_ str_
        where p !c !str 
                | c == '.' || c == ',' = parse ref (add cts ref) str
                | c == '^'             = parse ref (add cts ref) $ B.drop 1 str
                | c == '*' || c == '$' = parse ref cts str
                | c == '-' || c == '+' = let Just (cnt,rest) = B.readInt str
                                             var = (if c=='+' then Ins else Del) (B.unpack $ B.take (fromIntegral cnt) rest)
                                         in parse ref (addvar cts var) (B.drop (fromIntegral cnt) rest)
                | otherwise            = parse ref (add cts c) str
              
              add (C as cs gs ts vs) 'A' = (C (as+1) cs gs ts vs)
              add (C as cs gs ts vs) 'C' = (C as (cs+1) gs ts vs)              
              add (C as cs gs ts vs) 'G' = (C as cs (gs+1) ts vs)
              add (C as cs gs ts vs) 'T' = (C as cs gs (ts+1) vs)
              add (C as cs gs ts vs) 'N' = (C as cs gs ts vs)              
              add _ n = error ("Not a nucleotide: "++show n)
              addvar (C as cs gs ts vs) v = (C as cs gs ts (v:vs))

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
