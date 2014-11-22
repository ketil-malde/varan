-- Draw Unicode sparklines, using code points 0x2581 to 0x2588

module Main where

import System.Console.ANSI
import Data.Char (chr)
import Data.List (elemIndices, sort)
import Data.Maybe

import MPileup
import Count
import qualified Data.ByteString.Lazy.Char8 as BL

main = do
  inp <- BL.getContents
  let ms = map readPile1 $ BL.lines inp
  mapM_ putStrLn $ sparklines ms


test = putStrLn (concat [count2char [5,b,0,0] | b <- [0..10]]
                 ++concat [count2char [0,5,b,0] | b <- [0..10]]
                 ++ concat [count2char [0,0,5,b] |  b <- [0..10]] ++ setSGRCode [])

sparklines :: [MPileRecord] -> [String]
sparklines = map sparkline . transpose . map counts

transpose xs
  | any null xs = []
  | otherwise = map head xs : transpose (map tail xs)

sparkline :: [Counts] -> String
sparkline = (++setSGRCode []) . concatMap (count2char . toList)

-- display a frequency spectrum as a sparkline char
count2char [0,0,0,0] = setSGRCode [SetColor Background Dull Black] ++ " "
count2char counts = let
  ((a,b):(c,d):rest) = reverse $ sort $ zip counts [0..3]
  [(majpos,major),(minpos,minor)] = sort [(b,a),(d,c)]
  toCol 0 = Red -- A
  toCol 1 = Blue -- C
  toCol 2 = Green -- G 
  toCol 3 = Yellow -- T
  majfreq = 9 * major `div` (major + minor)
  in f2char (toCol majpos) (toCol minpos) majfreq
     -- show (majpos,minpos,majfreq,counts) -- 

-- From colors and frequency, calculate the appropriate sparkline char
-- f should range be in the range [0..9]
f2char :: Color -> Color -> Int -> String
f2char fg bg f 
  | f<=0 = cols ++ " " -- [chr 0x2581]
  | f>=9 = inv ++ " " -- [chr 0x2588]
  | otherwise = cols ++ [chr (0x2580+f)]
  where
    cols = setSGRCode [SetColor Foreground Dull fg,SetColor Background Dull bg]
    inv = setSGRCode [SetColor Foreground Dull bg,SetColor Background Dull fg]

reset :: IO ()
reset = setSGR []
