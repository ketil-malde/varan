{-# Language DeriveDataTypeable #-}

module Options (Options(..), getArgs) where

import System.Console.CmdArgs
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO

-- data SnpMode = Default | Pairwise | Groups deriving (Typeable,Data)

data Options = Opts 
  { suppress, variants
  , chi2, f_st, pi_k, conf, ds, esiv, pconf :: Bool
  , input, output :: FilePath
  , global :: Bool
  , threads :: Int
    , min_cov, max_cov :: Int
  } deriving (Typeable,Data)

defopts :: Options
defopts = Opts 
  { -- mode = Default &= help "mode of operation"
    output = "" &= help "output file name" &= typFile
  , suppress = False &= help "omit non-variant lines from output"
  , variants = False &= help "output list of non-SNP variants"
  , global = False &= help "calculate global statistics"
  , chi2   = False &= help "calculate chi² probability" &= ignore
  , f_st   = False &= help "estimate F_st" &= name "fst"
  , pi_k   = False &= help "estimate Pi_k"  
  , conf   = False &= help "check if major allele frequency confidence intervals overlap"
  , pconf  = False &= help "pairwise major allele confidence"
  , ds     = False &= help "output distance between major allele frequency confidence intervals"
  , esiv   = False &= help "output expected site information value for SNPs"
  , threads = 100 &= help "queue lenght for parallel execution"
  , min_cov = 0     &= help "minimum coverage to include"
  , max_cov = 32000 &= help "maximum coverage to include"
  , input  = [] &= args &= typFile
  } &= program "varan"
    &= summary "Identify genetic variants from pooled sequences."
    &= details ["Examples:",""
               ,"Read input from a pipe, calculate site-wise Fst and confidence intervals, ignoring non-variant sites:",""
               ,"  samtools mpileup -f ref.fasta | varan --fst --conf -s -o snps.txt",""
               ,"Read input from a file, send the site-wise output to /dev/null, and only output global statistics to standard output:",""
               ,"  varan --global -o /dev/null input.mpile",""
               ]    

getArgs :: IO (IO BL.ByteString,Options)
getArgs = do
  opts <- cmdArgs defopts
  inp  <- get_input opts
  return (inp,opts)
  
-- | Determine where and how to read the input mpile data.
get_input :: Options -> IO (IO BL.ByteString)
get_input o = do
  interactive <- hIsTerminalDevice stdin
  return $ if null (input o) then 
             if interactive then error "Cowardly refusing to read input from terminal.\nIf you really want this, specify '-' as input.\n(or use --help for help)."
             else BL.getContents
           else if input o == "-" then BL.getContents 
                else BL.readFile (input o)
