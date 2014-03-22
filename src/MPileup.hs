module MPileup (Counts(..), readPile1, toList, major_allele, by_major_allele, showC, showV, sumList, MPileRecord) where

import Data.Char (toUpper, isDigit)
import Data.List (foldl',intercalate,nub,elemIndex)
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
                              in (check trs, chr, pos, ref, map count trs)
    parse1 xs = error ("parse1: insufficiently long line:"++show xs)
    
    check ts = let t = concat ts in null t || all (==head t) t

    triples _ [] = []
    triples ref (_cnt:bases:_quals:rest) = parse ref (B.map toUpper bases) : triples ref rest
    triples _ _ = error "triples: incorrect number of columns"
    
    -- this could probalby be faster if counting was incorporated directly, 
    -- avoiding the intermediate [Variant] data structure
    parse :: Char -> B.ByteString -> [Variant]
    parse ref bs = case B.uncons bs of
      Nothing -> []
      Just (c_,str_) -> p c_ str_
        where p c str 
                | c == '^'             = parse ref $ B.drop 1 str
                | c == '*' || c == '$' = parse ref str                               
                | c == '.' || c == ',' = Nuc ref : parse ref str
                | c == '-' || c == '+' = let Just (cnt,rest) = B.readInt str
                                         in (if c=='+' then Ins else Del) (B.unpack $ B.take (fromIntegral cnt) rest) 
                                            : parse ref (B.drop (fromIntegral cnt) rest)
                | otherwise            = Nuc c : parse ref str

    count :: [Variant] -> Counts
    count vars = foldl' f (C 0 0 0 0 []) vars
      where f (C as cs gs ts vs) x = case x of
              Nuc 'A' -> (C (as+1) cs gs ts vs)
              Nuc 'C' -> (C as (cs+1) gs ts vs)
              Nuc 'G' -> (C as cs (gs+1) ts vs)
              Nuc 'T' -> (C as cs gs (ts+1) vs)
              Nuc 'N' -> (C as cs gs ts vs)
              Nuc _   -> error ("Not a nucleotide: "++show x++"\n"++show vars)
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
