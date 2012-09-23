module MPileup where

import Data.Char (toUpper, isDigit)

import AgrestiCoull

showPile :: (String,String,Char,[(Int,Int,Int,Int)]) -> String
showPile (chr,pos,ref,stats@(s1:ss)) = chr++"\t"++pos++"\t"++[ref]
                               ++concat ["\t"++show s | s <- stats]
                               ++"\t-"++concat ["\t"++conf s1 s | s <- ss] 

conf :: (Int, Int, Int, Int) -> (Int, Int, Int, Int) -> String
conf (a1,c1,g1,t1) (a2,c2,g2,t2) = let 
  s1 = a1+c1+g1+t1
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
readPile :: String -> [(String,String,Char,[(Int,Int,Int,Int)])]
readPile = map (parse1 . words) . lines
  where
    parse1 (chr:pos:(ref:_):rest) = (chr,pos,ref,triples ref rest)
    
    triples :: Char -> [String] -> [(Int,Int,Int,Int)]
    triples _ [] = []
    triples ref (_count:bases:_quals:rest) = count (filter (`elem` "ACGT") $ parse_bases ref bases) : triples ref rest
    
    parse_bases _ [] = []
    parse_bases ref (c:str) 
      | c == '^'             = parse_bases ref $ drop 1 str
      | c == '$'             = parse_bases ref str                               
      | c == '.' || c == ',' = ref : parse_bases ref str
      | c == '-' || c == '+' = let (cnt,rest) = span isDigit str
                               in parse_bases ref $ drop (read cnt) rest
      | otherwise            = toUpper c : parse_bases ref str
    count xs = (countCh 'A' xs, countCh 'C' xs, countCh 'G' xs, countCh 'T' xs) where countCh x = length . filter (==x)