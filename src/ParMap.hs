module ParMap (parMap) where

import Control.Concurrent
import Control.Monad
import System.IO.Unsafe

parMap :: Int -> (b -> a) -> [b] -> IO [a]
parMap t f xs = do
  outs <- replicateM t newEmptyMVar
  ins  <- replicateM t newEmptyMVar
  sequence_ [forkIO (worker i o) | (i,o) <- zip ins outs]
  forkIO $ sequence_ [putMVar i (f x) | (i,x) <- zip (cycle ins) xs]
  sequence' [takeMVar o | (o,_) <- zip (cycle outs) xs]

worker :: MVar a -> MVar a -> IO b
worker imv omv = forever $ do 
  ex <- takeMVar imv 
  ex `seq` putMVar omv ex

sequence' :: [IO a] -> IO [a]
sequence' [] = return []
sequence' (m:ms) = do
  x <- m
  xs <- unsafeInterleaveIO $ sequence' ms
  return (x:xs)