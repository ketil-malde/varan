module Main where

import MPileup (readPile1)
import Process (proc_fused, run_procs)
import Options
import Control.Monad (when)

import qualified Data.ByteString.Lazy.Char8 as BL

main :: IO ()
main = do
  (inp,o) <- Options.getArgs
  lns <- BL.lines `fmap` inp
  when (null lns) $ error "No lines in input!"
  if not (global o)
    then proc_fused o lns -- faster?
    else do
      let recs = map readPile1 lns -- seems faster for few CPUs?
      -- recs <- parMap readPile1 lns
      run_procs o recs
