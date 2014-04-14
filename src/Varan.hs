module Main where

import MPileup (readPile1)
import Process (proc_fused, run_procs, showPile')
import Options
import Control.Monad (when)
import ParMap

import qualified Data.ByteString.Lazy.Char8 as BL

main :: IO ()
main = do
  (inp,o) <- Options.getArgs
  lns <- BL.lines `fmap` inp
  when (null lns) $ error "No lines in input!"
  if not (global o)
    then proc_fused o lns -- faster?
    else do
      -- seems faster for few CPUs?
      let comb = showPile' o . readPile1 32000
      recs <- if threads o > 1 then parMap (threads o) comb lns 
              else return $ map comb lns
      run_procs o recs
