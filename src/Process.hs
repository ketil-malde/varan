module Process (proc_default, proc_fused, proc_globals) where

import Options
import ParMap
import MPileup
import Metrics

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import Text.Printf
import System.IO

type ProcF = Options -> [MPileRecord] -> IO ()

proc_default :: ProcF
proc_default o (l:ls) = do
  outh <- if null (output o) || output o == "-" then return stdout else openFile (output o) WriteMode
  B.hPutStr outh $ gen_header o l
  mapM_ (B.hPutStr outh) =<< parMap (threads o) (showPile o) (l:ls)
  hClose outh

proc_fused :: Options -> [BL.ByteString] -> IO ()
proc_fused o (l:ls) = do
  outh <- if null (output o) || output o == "-" then return stdout else openFile (output o) WriteMode
  B.hPutStr outh $ gen_header o $ readPile1 l
  mapM_ (B.hPutStr outh) =<< parMap (threads o) (showPile o . readPile1) (l:ls)
  hClose outh

proc_globals :: ProcF
proc_globals = undefined
               -- add to univar-like structure
               -- coverage stats, f_st, pi - and length.

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

