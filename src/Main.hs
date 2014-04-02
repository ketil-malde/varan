module Main where

import MPileup (readPile1)
import Process (proc_default)
import Options
import ParMap

import qualified Data.ByteString.Lazy.Char8 as BL

main :: IO ()
main = do
  (inp,o) <- Options.getArgs
  lns <- BL.lines `fmap` inp
  -- let inp = map readPile1 lns -- slow, serialized parsing
  recs <- parMap readPile1 lns
  case recs of 
    [] -> error "No lines in input!"
    (m:ms) -> proc_default o (m:ms)
