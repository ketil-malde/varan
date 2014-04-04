module Main where

import MPileup (readPile1)
import Process (proc_default, proc_fused)
import Options
import ParMap

import qualified Data.ByteString.Lazy.Char8 as BL

main :: IO ()
main = do
  (inp,o) <- Options.getArgs
  lns <- BL.lines `fmap` inp
  proc_fused o lns
{-  
  let recs = map readPile1 lns -- seems faster for few CPUs?
  -- recs <- parMap readPile1 lns
  case recs of 
    [] -> error "No lines in input!"
    (m:ms) -> proc_default o (m:ms)
-}
