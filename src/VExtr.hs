module Main where

import MPileup
import Count
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
        mp2char (MPR _ig _chr _cpos _ref cts) = let ss = ptSum cts
                                                    vs = map getV cts
                                                in selectChar ss
        selectChar ss = case toList ss of
          [0,0,0,0] -> 'n'
          [_,0,0,0] -> 'A'
          [0,_,0,0] -> 'C'
          [0,0,_,0] -> 'G'
          [0,0,0,_] -> 'T'

          [x,0,y,0] -> maybeLower x y 'a' 'g' 'R'
          [0,x,0,y] -> maybeLower x y 'c' 't' 'Y'
          [0,x,y,0] -> maybeLower x y 'c' 'g' 'S'
          [x,0,0,y] -> maybeLower x y 'a' 't' 'W'
          [0,0,x,y] -> maybeLower x y 'g' 't' 'K'
          [x,y,0,0] -> maybeLower x y 'a' 'c' 'M'

          [0,_,_,_] -> 'B'
          [_,0,_,_] -> 'D'
          [_,_,0,_] -> 'H'
          [_,_,_,0] -> 'V'

          _ -> 'N'
        maybeLower x y c1 c2 c3 = if x>1 && y >1 then c3 else if x>1 then c1 else c2
