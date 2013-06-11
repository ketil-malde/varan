module Main where

import System.Environment (getArgs)
import MPileup

main :: IO ()
main = do
  f <- getArgs
  ls <-  readPile `fmap` case f of 
    [] -> getContents 
    [fn] -> readFile fn
  mapM_ (\x -> do l <- showPile x
                  putStrLn l) ls
