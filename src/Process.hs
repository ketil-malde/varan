{-# Language BangPatterns #-}

module Process (proc_default, proc_fused, proc_globals, run_procs) where

import Options
import ParMap
import MPileup
import Metrics

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import Text.Printf
import System.IO
import Control.Concurrent

proc_fused :: Options -> [BL.ByteString] -> IO ()
proc_fused o (l:ls) = do
  outh <- if null (output o) || output o == "-" then return stdout else openFile (output o) WriteMode
  B.hPutStr outh $ gen_header o $ readPile1 l
  mapM_ (B.hPutStr outh) =<< parMap (threads o) (showPile o . readPile1) (l:ls)
  hClose outh

run_procs :: Options -> [MPileRecord] -> IO ()
run_procs o recs = do
  (li,lo) <- start_proc $ proc_default o
  (gi,go) <- start_proc $ proc_globals o
  let run (r:rs) = do
        push_procs (Just r) [li,gi]
        run rs
      run [] = do
        push_procs Nothing [li,gi]
        takeMVar lo -- collect is a no-op
        collect_globals =<< takeMVar go
  run recs

start_proc :: (MVar inv -> MVar out -> IO ()) -> IO (MVar inv, MVar out)
start_proc f = do
  imv <- newEmptyMVar
  omv <- newEmptyMVar
  forkIO $ f imv omv
  return (imv,omv)

push_procs :: a -> [MVar a] -> IO ()
push_procs r mvs = mapM_ (\m -> putMVar m r) mvs

proc_default :: Options -> MVar (Maybe MPileRecord) -> MVar () -> IO ()
proc_default o imv omv = do
  let use_stdout = null (output o) || output o == "-"
  outh <- if use_stdout then return stdout else openFile (output o) WriteMode
  Just l <- takeMVar imv -- or fail!
  B.hPutStr outh $ gen_header o l
  B.hPutStr outh $ showPile o l
  let run = do
        ml <- takeMVar imv
        case ml of
          Nothing -> do
            if (not use_stdout) then hClose outh else return ()
            putMVar omv ()
          Just x -> do
            B.hPutStr outh $ showPile o x
            run
  run

collect_globals :: Show a => a -> IO ()
collect_globals x = putStrLn ("Total nucleotide diversity: "++show x)

proc_globals :: Options -> MVar (Maybe MPileRecord) -> MVar Double -> IO ()
proc_globals _o imv omv = go 0.0
  where go !acc = do
          v <- takeMVar imv
          case v of
            Just (_,_,_,_,counts) -> do
              let p = Metrics.pi_k counts
              go (if isNaN p then acc else acc+p)
            Nothing -> putMVar omv acc

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

