module Main where

import MPileup (readPile1, showC, showV, by_major_allele, MPileRecord)
import Metrics(f_st, pi_k, conf_all, pearsons_chi², ds_all)
import Text.Printf
import Options
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO
import ParMap

main :: IO ()
main = do
  o <- Options.getArgs
  interactive <- hIsTerminalDevice stdin
  lns <- BL.lines `fmap` if null (input o) then 
                           if interactive then error "Cowardly refusing to read input from terminal.\nIf you really want this, specify '-' as input.\n(or use --help for help)."
                           else BL.getContents
                         else if input o == "-" then BL.getContents 
                              else BL.readFile (input o)
  outh <- if null (output o) || output o == "-" then return stdout else openFile (output o) WriteMode
  case lns of 
    [] -> error "No lines in input!"
    (l:ls) -> let hdr = gen_header o (readPile1 l)
              in do B.hPutStr outh hdr 
                    mapM_ (B.hPutStr outh) =<< parMap (showPile o . readPile1) (l:ls)
  hClose outh

-- generate the appropriate header, based on number of input pools
gen_header :: Options -> MPileRecord -> B.ByteString
gen_header o (_,_,_,_,cs) = B.pack $ concat [
  standard
  ,if Options.f_st o then "\tF_st" else ""
  ,if Options.pi_k o then "\tPi_k" else ""
  ,if Options.chi2 o then "\tChi²" else ""
  ,if Options.conf o then concat ["\tCI "++show n | n <- [1..(length cs)]] else ""
  ,if Options.conf o then "\tdelta-sigma" else ""
  ,if Options.variants o then "\tVariants" else ""
  ,"\n"  
  ]
  where
    standard = "#Target seq.\tpos\tref"++concat ["\tsample "++show x | x <- [1..(length cs)]]++"\tcover"

showPile :: Options -> MPileRecord -> B.ByteString
showPile _ (_,_,_,_,[]) = error "Pileup with no data?"
showPile o inp@(f,_,_,_,counts) = if suppress o && f then B.empty else (B.concat
          [ default_out inp
          , when (Options.f_st o) (printf "\t%.3f" (Metrics.f_st counts))
          , when (Options.pi_k o) (printf "\t%.3f" (Metrics.pi_k counts))
          , when (Options.chi2 o) (printf "\t%.3f" (Metrics.pearsons_chi² $ by_major_allele counts))
          , when (Options.conf o) (conf_all counts)
          , when (Options.ds o) ("\t"++(unwords $ map (\x -> if x>=0 then printf "%.2f" x else " -  ") $ ds_all 2.326 counts))
          , when (Options.variants o) ("\t"++showV counts)
          , B.pack "\n"
          ])
  where when p s = if p then B.pack s else B.empty
        
default_out :: MPileRecord -> B.ByteString
default_out (_,chr,pos,ref,stats) = 
  B.concat ([chr',tab,pos',tab,B.singleton ref]++samples++counts)
    where cnts = map showC stats
          tab  = B.pack "\t"
          samples = [B.append tab (B.pack s) | s <- map fst cnts]
          counts  = [tab,B.pack $ show $ sum $ map snd cnts] -- todo: add indels?
          chr' = B.concat (BL.toChunks chr)
          pos' = B.concat (BL.toChunks pos)

