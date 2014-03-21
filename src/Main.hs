module Main where

import MPileup (readPile1, showC, showV, by_major_allele, MPileRecord)
import Metrics(f_st, pi_k, conf_all, pearsons_chi²)
import Text.Printf
import Options
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO

main :: IO ()
main = do
  o <- Options.getArgs
  interactive <- hIsTerminalDevice stdin
  lns <- BL.lines `fmap` if null (input o) then 
                           if interactive then error "Cowardly refusing to read input from terminal.\nIf you really want this, specify '-' as input.\n(or use --help for help)."
                           else BL.getContents
                         else if input o == "-" then BL.getContents 
                              else BL.readFile (input o)
  let outf = if null (output o) || output o == "-" then BL.putStr else BL.writeFile (output o)
  case lns of 
    [] -> error "No lines in input!"
    (l:ls) -> let hdr = gen_header o (readPile1 l)
              in hdr `seq` mapM_ outf $ hdr : map (showPile o . readPile1) (l:ls)

-- generate the appropriate header, based on number of input pools
gen_header :: Options -> MPileRecord -> BL.ByteString
gen_header o (_,_,_,_,cs) = BL.pack $ concat [
  standard
  ,if Options.f_st o then "\tF_st" else ""
  ,if Options.pi_k o then "\tPi_k" else ""
  ,if Options.chi2 o then "\tChi²" else ""
  ,if Options.conf o then concat ["\tCI "++show n | n <- [1..(length cs)]] else ""
  ,if Options.variants o then "\tVariants" else ""
  ]
  where
    standard = "#Target seq.\tpos\tref"++concat ["\tsample "++show x | x <- [1..(length cs)]]++"\tcover"

showPile :: Options -> MPileRecord -> BL.ByteString
showPile _ (_,_,_,_,[]) = error "Pileup with no data?"
showPile o inp@(f,_,_,_,counts) = if suppress o && f then BL.empty 
                                  else BL.append (default_out inp) (BL.fromChunks 
          [ when (Options.f_st o) (printf "\t%.3f" (Metrics.f_st counts))
          , when (Options.pi_k o) (printf "\t%.3f" (Metrics.pi_k counts))
          , when (Options.chi2 o) (printf "\t%.3f" (Metrics.pearsons_chi² $ by_major_allele counts))
          , when (Options.conf o) (conf_all counts)
          , when (Options.variants o) ("\t"++showV counts)
          , B.pack "\n"
          ])
  where when p s = if p then B.pack s else B.empty
        
default_out :: MPileRecord -> BL.ByteString
default_out (_,chr,pos,ref,stats) = 
  BL.concat ([chr,tab,pos,tab,BL.singleton ref]++samples++counts)
    where cnts = map showC stats
          tab  = BL.pack "\t"
          samples = [BL.append tab (BL.pack s) | s <- map fst cnts]
          counts  = [tab,BL.pack $ show $ sum $ map snd cnts] -- todo: add indels?

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

