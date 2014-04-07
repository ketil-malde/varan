{-# Language BangPatterns #-}

module Process (proc_fused, run_procs) where

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
proc_fused _ [] = error "No input?"

-- | Runs a set of processes, distributes each MPileRecord to them
--   and runs the finalizer (collecting and outputting the results)  
run_procs :: Options -> [MPileRecord] -> IO ()
run_procs o recs = do
  (li,lfin) <- start_proc (proc_default o) return
  (gi,gfin) <- start_proc proc_gpi
               (\x -> putStrLn ("Global pi_k (nucleotide diversity): "++show x++"\n"))
  (ppi,pfin) <- start_proc proc_gppi out_gppi
  (fi,ffin) <- start_proc proc_gfst out_gfst
  let run (r:rs) = do
        push_procs (Just r) [li,gi,ppi,fi]
        run rs
      run [] = do
        push_procs Nothing [li,gi,ppi,fi]
        sequence_ [lfin,putStrLn "",gfin,pfin,ffin]
  run recs

-- | The main process (in the first parameter) reads from 'inv' and puts the result
--   in 'out' when it's done. It returns the input MVar, and the IO action that
--   runs the finalizer (second parameter) on the final result. 
start_proc :: (MVar inv -> MVar out -> IO ()) -> (out -> IO ()) -> IO (MVar inv, IO ())
start_proc f g = do
  imv <- newEmptyMVar
  omv <- newEmptyMVar
  _ <- forkIO $ f imv omv
  return (imv, takeMVar omv >>= g)

-- | Shortcut for feeding a datum to a list of process inputs.
push_procs :: a -> [MVar a] -> IO ()
push_procs r mvs = mapM_ (\m -> putMVar m r) mvs

proc_fold :: a -> (MPileRecord -> a -> a) -> MVar (Maybe MPileRecord) -> MVar a -> IO ()
proc_fold z f mv e = go z 
  where go !c = do
          v <- takeMVar mv
          case v of
            Nothing -> putMVar e c 
            Just x  -> go (f x c)

-- --------------------------------------------------
-- Actual output calculations

-- | Calculate global nucleotide diversity
proc_gpi :: MVar (Maybe MPileRecord) -> MVar Double -> IO ()
proc_gpi = proc_fold 0.0 f
  where f mpr cur =
          let p = Metrics.pi_k (counts mpr)
          in if ignore mpr || isNaN p then cur else cur+p

-- | Collect variation within and between for global pairwise Fst
proc_gfst :: MVar (Maybe MPileRecord) -> MVar [[(Double, Double)]] -> IO ()
proc_gfst = proc_fold zero f
  where f (MPR sup _ _ _ cts) cur =
          let new = Metrics.fst_params cts
          in if sup then cur else deepSeq $ zipWith (zipWith plus) cur new
        zero = repeat (repeat (0,0))
        plus (a,c) (b,d) = (a+b,c+d)
        deepSeq x | x == x = x

-- | Calculate and print global pairwise Fst
out_gfst :: [[(Double,Double)]] -> IO ()
out_gfst xs = do
  putStrLn "Pairwise FST values:"
  putStrLn (" "++ concat [ "    s"++show (i+1) | i <- [1..length $ head xs]])
  go 1 xs
  where go i (l:ls) = do 
          putStr ("s"++show i++replicate (6*i-4) ' ')
          putStrLn $ unwords $ map (\(t,w) -> printf "%.3f" ((t-w)/t)) l
          go (i+1) ls
        go _ [] = return ()

-- --------------------------------------------------

data UniVar = UV { count :: !Int, sum1, sum2 :: !Double }

add_uv :: UniVar -> Double -> UniVar
add_uv (UV c s s2) d = UV (c+1) (s+d) (s2+d*d)

-- | Collect nucleotide diversity within and between for global pairwise ND (pi)
proc_gppi :: MVar (Maybe MPileRecord) -> MVar (UniVar,[[Double]]) -> IO ()
proc_gppi = proc_fold zero f
  where f (MPR sup _ _ _ cts) (uv,cur) =
          let new = Metrics.ppi_params cts
              cov = fromIntegral $ sum $ map (\(C a c g t vs) -> a+c+g+t+length vs) cts
              nu = add_uv uv cov
              nc = if sup then cur else deepSeq $ zipWith (zipWith plus) cur new 
          in nu `seq` nc `seq` (nu,nc )
        zero = (UV 0 0 0, repeat (repeat 0))
        plus a b = if isNaN b then a else a+b
        deepSeq x | x == x = x


-- | Calculate and print global pairwise ND
--   Todo: divide by genome size        
out_gppi :: (UniVar,[[Double]]) -> IO ()
out_gppi (UV n s1 s2,xs) = do
  let n' = fromIntegral n
      go i (l:ls) = do 
          let s = show i in putStr ("s"++s++replicate (7*i-3-length s) ' ')
          putStrLn $ unwords $ map (\t -> printf "%.4f" (t/n')) l
          go (i+1) ls
      go _ [] = return ()
  putStrLn "Coverage statistics:"
  printf "  covered sites: %d\n  avg. cover: %.2f\n  std. dev.: %.2f\n\n" n  (s1/n') (sqrt (s2/n'-(s1*s1)/(n'*n')))
  putStrLn "Pairwise Nucleotide Diversities:"
  putStrLn (" "++ concat [ "     s"++show i | i <- [1..length (head xs)]])
  go 1 xs
  putStrLn ""

-- | Calculating per site statistics writing to a file or standard out.
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

-- generate the appropriate header, based on number of input pools
gen_header :: Options -> MPileRecord -> B.ByteString
gen_header o (MPR _ _ _ _ cs) = B.pack $ concat [
  standard
  ,if Options.f_st o then "\tF_st" else ""
  ,if Options.pi_k o then "\tPi_k" else ""
  ,if Options.chi2 o then "\tChi²" else ""
  ,if Options.conf o then concat ["\tCI "++show n | n <- [1..(length cs)]] else ""
  ,if Options.ds o then "\tdelta-sigma" else ""
  ,if Options.variants o then "\tVariants" else ""
  ,"\n"  
  ]
  where
    standard = "#Target seq.\tpos\tref"++concat ["\tsample "++show x | x <- [1..(length cs)]]++"\tcover"

showPile :: Options -> MPileRecord -> B.ByteString
showPile _ (MPR _ _ _ _ []) = error "Pileup with no data?"
showPile o mpr = if suppress o && ignore mpr then B.empty else (B.concat
          [ default_out mpr
          , when (Options.f_st o) (printf "\t%.3f" (Metrics.f_st $ counts mpr))
          , when (Options.pi_k o) (printf "\t%.3f" (Metrics.pi_k $ counts mpr))
          , when (Options.chi2 o) (printf "\t%.3f" (Metrics.pearsons_chi² $ by_major_allele $ counts mpr))
          , when (Options.conf o) (conf_all $ counts mpr)
          , when (Options.ds o) ("\t"++(unwords $ map (\x -> if x>=0 then printf "%.2f" x else " -  ") $ ds_all 2.326 $ counts mpr))
          , when (Options.variants o) ("\t"++showV (counts mpr))
          , B.pack "\n"
          ])
  where when p s = if p then B.pack s else B.empty
        
-- | The default output, with only coverage statistics
default_out :: MPileRecord -> B.ByteString
default_out (MPR _ chr pos ref stats) =
  B.concat ([chr',tab,pos',tab,B.singleton ref]++samples++fmtcounts)
    where cnts = map showC stats
          tab  = B.pack "\t"
          samples = [B.append tab (B.pack s) | s <- map fst cnts]
          fmtcounts  = [tab,B.pack $ show $ sum $ map snd cnts] -- todo: add indels?
          chr' = B.concat (BL.toChunks chr)
          pos' = B.concat (BL.toChunks pos)

