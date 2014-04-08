module Main where

import MPileup (readPile1)
import Process (proc_fused, run_procs)
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
      recs <- if threads o > 1 then parMap (threads o) readPile1 lns 
              else return $ map readPile1 lns
      run_procs o recs
