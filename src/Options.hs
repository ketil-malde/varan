{-# Language DeriveDataTypeable #-}

module Options (Options(..), getArgs) where

import System.Console.CmdArgs
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO

-- data SnpMode = Default | Pairwise | Groups deriving (Typeable,Data)

data Options = Opts 
  { suppress, variants
  , chi2, f_st, pi_k, conf, ds :: Bool
  , input, output :: FilePath
  } deriving (Typeable,Data)

defopts :: Options
defopts = Opts 
  { -- mode = Default &= help "mode of operation"
    output = "" &= help "output file name" &= typFile
  , suppress = False &= help "omit non-variant lines from output"
  , variants = False &= help "output list of non-SNP variants"
  , chi2   = False &= help "calculate chi² probability"
  , f_st   = False &= help "estimate F_st" &= name "fst"
  , pi_k   = False &= help "estimate Pi_k"  
  , conf   = False &= help "check if major allele frequency confidence intervals overlap"
  , ds     = False &= help "output distance between major allele frequency confidence intervals"
  , input  = [] &= args &= typFile
  } &= details ["Identify variants from mpileup data"]

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
