{-# Language GADTs #-}

module Variants where
import Data.Char (isDigit)

data Counts = C !Int !Int !Int !Int [Variant]

-- Convert 'Counts' to a list of allele counts
toList :: Num a => Counts -> [a]
toList (C x y z v _) = map fromIntegral [x,y,z,v]

-- | Pointwise summation of the input lists
sumList :: Num a => [[a]] -> [a]
sumList = foldr (zipWith (+)) [0,0,0,0]


data Variant where
  Nuc :: Char -> Variant
  Ins, Del :: String -> Variant 
  deriving (Show,Eq)

-- manually inlined in the `MPileup` module.  Don't use this one.
parse :: Char -> String -> [Variant]
parse _ [] = []
parse ref (c:str) 
      | c == '^'             = parse ref $ drop 1 str
      | c == '*' || c == '$' = parse ref str                               
      | c == '.' || c == ',' = Nuc ref : parse ref str
      | c == '-' || c == '+' = let (x,rest) = span isDigit str
                                   cnt = read x
                               in (if c=='+' then Ins else Del) (take cnt rest) 
                                  : parse ref (drop cnt rest)
      | otherwise            = Nuc c : parse ref str

