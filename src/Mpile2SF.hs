{-
  Convert MPile to SweepFinder input, i.e.
  <position> <x> <n> <folded>

  I guess all will have to be folded?  Or use ref as ancestor?  How does it matter?
-}

module Main where

import qualified Data.ByteString.Lazy.Char8 as B
import Text.Printf
import Data.Maybe

import MPileup
import Count

main :: IO ()
main = do
  mps <- map readPile1 `fmap` B.lines `fmap` B.getContents
  putStrLn "position\tx\tn\tfolded"
  mapM_ B.putStrLn . catMaybes $ map mp2sf mps

mp2sf :: MPileRecord -> Maybe B.ByteString
mp2sf mpr
  | ignore mpr = Nothing
  | not (refnuc mpr `elem` "ACGTacgt") = Nothing
  | otherwise  = Just (B.append (cpos mpr) (B.pack (printf "\t%d\t%d\t1"  refs tot)))
  where alleles = toList $ ptSum $ counts mpr :: [Int]
        tot  = sum alleles :: Int
        refs = alleles!!(case refnuc mpr of {'A' -> 0; 'C' -> 1; 'G' -> 2; 'T' -> 3; 'a' -> 0; 'c' -> 1; 'g' -> 2; 't' -> 3}) :: Int
