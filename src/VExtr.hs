module Main where

import MPileup
import Count
import Data.List (foldl')
import qualified Data.ByteString.Lazy.Char8 as BL

-- TODO: merge as separate mode

main :: IO ()
main = do
  inp <- BL.getContents
  let ms = map readPile1 $ BL.lines inp
  putStr (makeConsensus ms)

makeConsensus :: [MPileRecord] -> String
makeConsensus = map mp2char
  where mp2char :: MPileRecord -> Char
        mp2char (MPR ig chr cpos ref cts) = let ss = ptSum cts
                                                vs = map getV cts
                                            in selectChar ss
        selectChar ss = case toList ss of
          [0,0,0,0] -> 'n'
          [x,0,0,0] -> 'A'
          [0,x,0,0] -> 'C'
          [0,0,x,0] -> 'G'
          [0,0,0,x] -> 'T'

          [x,0,y,0] -> 'R'
          [0,x,0,y] -> 'Y'
          [0,x,y,0] -> 'S'
          [x,0,0,y] -> 'W'
          [0,0,x,y] -> 'K'
          [x,y,0,0] -> 'M'
          [0,x,y,z] -> 'B'
          [x,0,y,z] -> 'D'
          [x,y,0,z] -> 'H'
          [x,y,z,0] -> 'V'

          _ -> 'N'
