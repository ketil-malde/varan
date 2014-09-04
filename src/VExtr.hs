{-# Language DeriveDataTypeable #-}

module Main where

import MPileup
import Count
import qualified Data.ByteString.Lazy.Char8 as BL

import System.Console.CmdArgs

data Options = Opts { infile, outfile :: Maybe FilePath
                    , xout :: Bool
                    , mincount :: Int } deriving (Data,Typeable)

opts :: Options
opts = Opts 
  { infile = Nothing &= args &= typFile
  , outfile = Nothing &= help "output file"
  , xout = False &= help "output X instead of IUPAC codes for variable sites"
  , mincount = 1 &= help "ignore counts less than this"
  } &= program "vextr v0.4"
    &= summary "Extract consensus sequence from pooled sequences"

main :: IO ()
main = do
  o <- cmdArgs opts
  inp <- case infile o of Nothing -> BL.getContents
                          Just f -> BL.readFile f
  let ms = map readPile1 $ BL.lines inp
  (case outfile o of Nothing -> putStr
                     Just f -> writeFile f) (makeConsensus (xout o,mincount o) ms)

makeConsensus :: (Bool,Int) -> [MPileRecord] -> String
makeConsensus (iup,mct) = if iup then map (fixiup . mp2char) else map mp2char
  where maybeLower x y c1 c2 c3 = if x>mct && y >mct then c3 else if x>mct then c1 else c2
        fixiup c | c `elem` "ACGTacgtNn" = c
                 | otherwise             = 'X'
        mp2char :: MPileRecord -> Char
        mp2char (MPR _ig _chr _cpos _ref cts) = let ss = ptSum cts
                                                    -- vs = map getV cts
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
