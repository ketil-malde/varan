module Main where

import System.Environment (getArgs)
import MPileup

main :: IO ()
main = do
  f <- getArgs
  ls <-  readPile `fmap` case f of 
    [] -> getContents 
    [fn] -> readFile fn
  putStrLn $ unlines $ map showPile ls
