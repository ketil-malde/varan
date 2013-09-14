module Main where

import System.Environment (getArgs)
import MPileup (Counts, readPile, showC, showV, major_allele) 
import Metrics(angle, f_st, pi_k, conf, ci_dist, delta_sigma, pearson_chi²)
import Text.Printf
{-
import RandomSelect
import System.Random
import Control.Parallel
-}

main :: IO ()
main = do
  f <- getArgs
  ls <-  readPile `fmap` case f of 
    [] -> getContents 
    [fn] -> readFile fn
  gen_header $ dataids (head ls)
  mapM_ (\x -> do l <- showPile x
                  putStrLn l) ls

-- generate the appropriate header, based on number of input pools
gen_header :: [String] -> IO ()
gen_header _ = return ()

dataids :: (String,String,Char,[Counts]) -> [String]
dataids _ = []

showPile :: (String,String,Char,[Counts]) -> IO String
showPile (_,_,_,[]) = error "Pileup with no data?"
showPile (chr,pos,ref,stats@(s1:ss)) = do
{-  g <- newStdGen
  let (f,pf) = pval g f_st (s1:ss)
      (p,pp) = pval g pi_k (s1:ss)
  pf `par` pp `pseq` -}
  let cnts = map showC stats
  return (
    chr++"\t"++pos++"\t"++[ref]++concat ["\t"++s | s <- map fst cnts]++"\t"++show (sum $ map snd cnts)
    ++concat ["\t"++conf s1 s | s <- ss]
    --  ++ conf_all (s1:ss)
    ++concat [printf "\t%.3f" (angle s1 s) | s <- ss]
    ++printf "\t%.3f" (f_st (s1:ss))  -- print_pval (f,pf)
    ++printf "\t%.3f" (pi_k (s1:ss))  -- print_pval (p,pp)
    ++concat [printf "\t%.2f" (uncurry ci_dist $ major_allele s1 s) | s <- ss]
    ++let as = [major_allele s1 s | s <- ss] in printf "\t%.2f" (pearson_chi² (fst (head as):map snd as))
    ++concat [printf "\t%.2f" (uncurry (delta_sigma 1) $ major_allele s1 s) | s <- ss]
    ++concat [printf "\t%.2f" (uncurry (delta_sigma 2) $ major_allele s1 s) | s <- ss]
    ++"\t"++showV stats)

{-
print_pval :: (Double, Double) -> String
print_pval (a,b) = printf "\t%.3f p=%.3f" a b
-}

