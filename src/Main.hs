module Main where

import MPileup (Counts, readPile, showC, showV, by_major_allele)
import Metrics(f_st, pi_k, conf_all, ci_dist, delta_sigma, pearsons_chi²)
import Text.Printf
import Options

{-
import RandomSelect
import System.Random
-}

main :: IO ()
main = do
  o <- Options.getArgs
  ls <-  readPile `fmap` if null (input o) then getContents else readFile (input o)
  gen_header o (head ls)
  mapM_ (putStr . showPile o) ls

-- generate the appropriate header, based on number of input pools
gen_header :: Options -> (Bool,String,String,Char,[Counts]) -> IO ()
gen_header o (_,_,_,_,cs) = putStrLn $ concat [
  standard
  ,if Options.f_st o then "\tF_st" else ""
  ,if Options.pi_k o then "\tPi_k" else ""
  ,if Options.chi2 o then "\tChi²" else ""
  ,if Options.conf o then "\tAG-CI" else ""
  ]
  where
    standard = "#Target seq.\tpos\tref"++concat ["\tsample "++show x | x <- [1..(length cs)]]++"\tcover"

showPile :: Options -> (Bool,String,String,Char,[Counts]) -> String
showPile _ (_,_,_,_,[]) = error "Pileup with no data?"
showPile o inp@(f,_,_,_,counts) = if suppress o && f then "" else concat [
          default_out inp
          , if Options.f_st o then printf "\t%.3f" (Metrics.f_st counts) else ""
          , if Options.pi_k o then printf "\t%.3f" (Metrics.pi_k counts) else ""
          , if Options.chi2 o then printf "\t%.3f" (Metrics.pearsons_chi² $ by_major_allele counts) else ""
          , if Options.conf o then conf_all counts else ""  
          -- , if Options.ds o   then ?
          ,"\n" 
          ]

default_out :: (Bool,String, String, Char, [Counts]) -> String
default_out (_,chr,pos,ref,stats) = 
  chr++"\t"++pos++"\t"++[ref]++concat ["\t"++s | s <- map fst cnts]++"\t"++show (sum $ map snd cnts) -- todo: add indels?
  where cnts = map showC stats

{-  Confidence intervals?   
    ++concat ["\t"++conf s1 s | s <- ss]
    --  ++ conf_all (s1:ss)
-}
    
{- angle - not so useful, I think?
  ++concat [printf "\t%.3f" (angle s1 s) | s <- ss]
-}

{-
    ++concat [printf "\t%.2f" (uncurry ci_dist $ major_allele s1 s) | s <- ss]
-}
    
{-  ++concat [printf "\t%.2f" (uncurry (delta_sigma 1) $ major_allele s1 s) | s <- ss]
    ++concat [printf "\t%.2f" (uncurry (delta_sigma 2) $ major_allele s1 s) | s <- ss]
    ++"\t"++showV stats
-}

{-
print_pval :: (Double, Double) -> String
print_pval (a,b) = printf "\t%.3f p=%.3f" a b
-}

