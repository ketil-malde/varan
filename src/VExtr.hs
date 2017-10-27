{-# Language DeriveDataTypeable #-}

module VExtr where

import MPileup
import Count
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Char (toLower,toUpper)

import System.Console.CmdArgs
import Options (version, citation)

data Options = Opts { infile, outfile :: Maybe FilePath
                    , format :: Format
                    , fasta, pad :: Bool
                    , mincount :: Int
                    , minfreq :: Int } deriving (Data,Typeable)

opts :: Options
opts = Opts 
  { infile = Nothing &= args &= typFile
  , outfile = Nothing &= help "output file"
  , format = IUPAC &= help "output X, N, or [a/b] instead of IUPAC codes for variable sites"
  , fasta  = False &= help "output FASTA header"
  , pad    = False &= help "pad start of output with Ns"
  , mincount = 1 &= help "ignore counts less than this"
  , minfreq  = 5 &= help "ignore allele frequencies less than this (in percent)"
  } &= program ("vextr "++version)
    &= summary "Extract consensus sequence from pooled sequences"
    &= details (["Examples:", ""
             ,"Read input from a pipe, output IUPAC codes:"
             ,"", "  samtools mpileup -f ref.fasta reads.bam | vextr --format=iupac", ""
             ,"Read input from a file, create consensus FASTA sequence:"
             ,"", "  vextr input.mpile --fasta -o output.fasta", ""
             ]++citation)

data Format = Xs | IUPAC | Regex deriving (Data,Typeable,Show)

main :: IO ()
main = do
  o <- cmdArgs opts
  inp <- case infile o of Nothing -> BL.getContents
                          Just f -> BL.readFile f
  let ms = map readPile1 $ BL.lines inp
      outf = case outfile o of Nothing -> putStr
                               Just f -> writeFile f
      gen = if fasta o then makeFasta else makeConsensus
  outf $ gen (format o,mincount o,minfreq o,pad o) ms


makeFasta :: (Format,Int,Int,Bool) -> [MPileRecord] -> String
makeFasta fi ms = let
  header = case ms of
    (m1:_) -> '>':BL.unpack (chrom m1)++posinf
      where posinf = case fi of
              (_,_,_,False) -> ":"++BL.unpack (cpos m1)
              _ -> ""
    [] -> ""
  breaks str = case splitAt 60 str of
    (rest,"") -> [rest]
    (this,more) -> this : breaks more
  in unlines (header:breaks (makeConsensus fi ms))

makeConsensus :: (Format,Int,Int,Bool) -> [MPileRecord] -> String
makeConsensus (iup,mct,mfq,pd) ms = prefix ++ concatMap (fixiup iup . selectChar mct mfq . ptSum . counts) ms
  where prefix = if pd then replicate padlen 'N' else []
        padlen = case ms of
          (m1:_) -> read (BL.unpack (cpos m1)) - 1
          _ -> 0

-- | Optionally change from IUPAC code to X or regex
fixiup :: Format -> Char -> String
fixiup iup c | c `elem` "ACGTacgt" = [c]
             | otherwise             = case iup of
          Xs -> "X"
          IUPAC -> [c]
          -- this should use upper/lower case to show (very) minor alleles
          Regex -> case toUpper c of
            'R' -> "[A/G]"
            'Y' -> "[C/T]"
            'S' -> "[C/G]"
            'W' -> "[A/T]"
            'K' -> "[G/T]"
            'M' -> "[A/C]"
            'B' -> "[C/G/T]"
            'D' -> "[A/G/T]"
            'H' -> "[A/C/T]"
            'V' -> "[A/C/G]"
            'N' -> "[A/C/G/T]"
            x   -> [x]

-- | Convert allele counts into IUPAC character
selectChar :: Int -> Int -> Counts -> Char
selectChar mct mfq ss =
  let xs = toList ss
      t  = sum xs
      as = map (>0) xs
      bs = map (\x -> x >= mct && x >= t*mfq`div`100) xs
      char = case bs of
          [False,False,False,False] -> 'n'
          [True ,False,False,False] -> 'A'
          [False,True, False,False] -> 'C'
          [False,False,True, False] -> 'G'
          [False,False,False,True ] -> 'T'

          [True, False,True ,False] -> 'R'
          [False,True ,False,True ] -> 'Y'
          [False,True ,True ,False] -> 'S'
          [True ,False,False,True ] -> 'W'
          [False,False,True ,True ] -> 'K'
          [True ,True ,False,False] -> 'M'

          [False,True ,True ,True ] -> 'B'
          [True ,False,True ,True ] -> 'D'
          [True ,True ,False,True ] -> 'H'
          [True ,True ,True ,False] -> 'V'
          _ -> 'N'
      in (if as == bs then id else toLower) char

{-
   a better way to do this, might be to have selectChar return the set of alleles
   in upper/lower case, sorted by frequency.  So e.g. "GTa" means that G and T are present
   above some threshold, a is rare, and c is seen rarely enough to be ignored.  IUPAC might then be 'k',
   and regex might be [G/T/a].  Do other programs expect regexps to be ordered?

   Another thing to consider might be multiple samples.  A couple of rare alleles in a single sample
   are likely read errors, but if they occur in multiple samples, they are more likely to be real.
   OTOH, errors depend on local context, so errors in samples are not independent... not so easy, this.
-}
