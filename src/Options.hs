{-# Language DeriveDataTypeable #-}

module Options (Options(..), getArgs) where

import System.Console.CmdArgs

-- data SnpMode = Default | Pairwise | Groups deriving (Typeable,Data)

data Options = Opts 
  { suppress
  , chi2, f_st, pi_k, conf, ds :: Bool
  , input, output :: FilePath
  } deriving (Typeable,Data)

defopts :: Options
defopts = Opts 
  { -- mode = Default &= help "mode of operation"
    output = "" &= help "output file name" &= typFile
  , suppress = False &= help "omit non-variant lines from output"
  , chi2   = False &= help "calculate chi² probability"
  , f_st   = False &= help "estimate F_st"
  , pi_k   = False &= help "estimate Pi_k"  
  , conf   = False &= help "check if major allele frequency confidence intervals overlap"
  , ds     = False &= help "output distance between major allele frequency confidence intervals"
  , input  = [] &= args &= typFile
  } &= details ["Identify variants from mpileup data"]

getArgs :: IO Options
getArgs = cmdArgs defopts
