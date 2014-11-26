{-# Language DeriveDataTypeable #-}

module VExtr where

import MPileup
import Count
import qualified Data.ByteString.Lazy.Char8 as BL

import System.Console.CmdArgs

data Options = Opts { infile, outfile :: Maybe FilePath
                    , format :: Format
                    , fasta :: Bool
                    , mincount :: Int
                    , minfreq :: Int } deriving (Data,Typeable)

opts :: Options
opts = Opts 
  { infile = Nothing &= args &= typFile
  , outfile = Nothing &= help "output file"
  , format = IUPAC &= help "output X, N, or [a/b] instead of IUPAC codes for variable sites"
  , fasta  = False &= help "output FASTA header"
  , mincount = 1 &= help "ignore counts less than this"
  , minfreq  = 5 &= help "ignore allele frequencies less than this"
  } &= program "vextr v0.5"
    &= summary "Extract consensus sequence from pooled sequences"
    &= details ["Examples:", ""
             ,"Read input from a pipe, output IUPAC codes:"
             ,"", "  samtools mpileup -f ref.fasta reads.bam | vextr --format=iupac", ""
             ,"Read input from a file, create consensus FASTA sequence:"
             ,"", "  vextr input.mpile --fasta -o output.fasta", ""
             ,"If you use this program, please cite:"
             ,"  BMC Genomics 2014, 15(Suppl 6):S20"
             ,"  http://www.biomedcentral.com/1471-2164/15/S6/S20"
               ]

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
  outf $ gen (format o,mincount o,minfreq o) ms


makeFasta :: (Format,Int,Int) -> [MPileRecord] -> String
makeFasta fi ms = let
  header = case ms of
    (m1:_) -> BL.unpack (chrom m1)++":"++BL.unpack (cpos m1)
    [] -> ""
  breaks str = case splitAt 60 str of
    (rest,"") -> [rest]
    (this,more) -> this : breaks more
  in unlines (header:breaks (makeConsensus fi ms))

makeConsensus :: (Format,Int,Int) -> [MPileRecord] -> String
makeConsensus (iup,mct,mfq) = concatMap (fixiup iup . selectChar mct mfq . ptSum . counts)
  -- todo: include variants

-- this doesn't work so well with high coverage/many libraries

-- | Optionally change from IUPAC code to X or regex
fixiup :: Format -> Char -> String
fixiup iup c | c `elem` "ACGTacgtNn" = [c]
             | otherwise             = case iup of
          Xs -> "X"
          IUPAC -> [c]
          Regex -> case c of
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
            x   -> [x]

-- | Convert allele counts into IUPAC character
selectChar :: Int -> Int -> Counts -> Char
selectChar mct mfq ss = case toList ss of
          [0,0,0,0] -> 'n'
          [_,0,0,0] -> 'A'
          [0,_,0,0] -> 'C'
          [0,0,_,0] -> 'G'
          [0,0,0,_] -> 'T'

          [x,0,y,0] -> maybeWild x y 'a' 'g' 'R'
          [0,x,0,y] -> maybeWild x y 'c' 't' 'Y'
          [0,x,y,0] -> maybeWild x y 'c' 'g' 'S'
          [x,0,0,y] -> maybeWild x y 'a' 't' 'W'
          [0,0,x,y] -> maybeWild x y 'g' 't' 'K'
          [x,y,0,0] -> maybeWild x y 'a' 'c' 'M'

          [0,_,_,_] -> 'B'
          [_,0,_,_] -> 'D'
          [_,_,0,_] -> 'H'
          [_,_,_,0] -> 'V'

          _ -> 'N'
  where maybeWild x y c1 c2 c3 =
          let xok = x > mct && x > (x+y)*mfq`div`100
              yok = y > mct && y > (x+y)*mfq`div`100
          in if xok && yok then c3 else if xok then c1 else if yok then c2 else 'n'
