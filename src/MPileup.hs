{-# Language BangPatterns #-}
module MPileup (Counts(..), readPile1, toList, major_allele, by_major_allele, showC, showV, sumList, MPileRecord(..)) where

import Data.Char (toUpper)
import Data.List (intercalate,nub,elemIndex)
import qualified Data.ByteString.Lazy.Char8 as B
import Variants hiding (parse)
import Count

data MPileRecord = MPR { ignore :: !Bool
                       , chrom, cpos :: !B.ByteString
                       , refnuc :: !Char
                       , counts :: ![Counts]
                       }

-- convert counts to major/non-major allele counts
by_major_allele :: [Counts] -> [[Int]] -- always length 2
by_major_allele cs = let
  ls = map (toList . getcounts) cs
  s  = sumList ls
  Just i = elemIndex (maximum s) s
  in map (\l -> [l!!i,sum l-l!!i]) ls

-- pick out major allele in first count, and output number of same/different
major_allele :: Counts -> Counts -> ((Int,Int),(Int,Int))
major_allele (C x _) (C y _) =
  let s1 = toList x
      s2 = toList y
      m = maximum s1
      Just i = elemIndex m s1
  in ((m, sum s1-m),(s2!!i,sum [s2!!j | j <- [0..3], j /= i]))

-- count major allele in first sample
-- return flag whether informative, chrom, pos, ref, 
readPile1 :: B.ByteString -> MPileRecord 
readPile1 = parse1 . B.split '\t'  -- later samtools sometimes outputs empty strings in columns
  where
    parse1 (chr:pos:r:rest) = let trs = triples ref rest
                                  ref = B.head r
                              in MPR (ign trs) chr pos ref trs
    parse1 xs = error ("parse1: insufficiently long line:"++show xs)
    
    ign cs = (<=1) $ length $ filter (/=(0::Int)) $ toList $ ptSum $ map getcounts cs

    triples _ [] = []
    triples ref (_cnt:bases:_quals:rest) = let this = parse ref (C 0 []) (B.map toUpper bases) 
                                                      in this `seq` this : triples ref rest
    triples _ _ = error "triples: incorrect number of columns"
    
    -- note: does not (yet) deal with '<' and '>', which apparently can occur
    -- perhaps it would be faster to do this as a sequence of BS ops, which would fuse?
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
              
              add (C x vs) 'A' = (C (addA x 1) vs)
              add (C x vs) 'C' = (C (addC x 1) vs)
              add (C x vs) 'G' = (C (addG x 1) vs)
              add (C x vs) 'T' = (C (addT x 1) vs)
              add c 'N' = c
              add _ n = error ("Not a nucleotide: "++show n)
              addvar (C x vs) v = (C x (v:vs))

-- | Show SNP counts and coverage
showC :: Counts -> (String,Int)
showC (C x _) = (" "++(intercalate ":" $ map show (toList x :: [Int])),covC x)

-- | Show structural variant count
showV :: [Counts] -> String
showV cs = let
  getv (C _ v) = v
  vs = nub $ concatMap getv cs
  countV :: Variant -> Counts -> Int
  countV v c = length . filter (==v) $ getv c
  in intercalate "\t" (show vs:[unwords $ map (\v -> show $ countV v c) vs | c <- cs])
