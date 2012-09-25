{-# Language GADTs #-}

module Variants where
import Data.Char (isDigit)

data Variant where
  Nuc :: Char -> Variant
  Ins, Del :: String -> Variant 
  deriving (Show,Eq)

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

