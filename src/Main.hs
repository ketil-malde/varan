module Main where

import System.Environment (getArgs)
import MPileup

main :: IO ()
main = do
  [f] <- getArgs
  ls <-  readPile `fmap` readFile f
  putStrLn $ unlines $ map showPile ls
