{-# Language DeriveDataTypeable #-}

-- Draw Unicode sparklines, using code points 0x2581 to 0x2588

module Main where

import System.Console.ANSI
import Data.Char (chr)
import Data.List (sort)
import System.Console.CmdArgs

import ESIV
import Count
import VExtr (makeConsensus, Format(IUPAC))
import MPileup
-- import Count
import qualified Data.ByteString.Lazy.Char8 as BL

data Options = Test
             | Disp
             | Info -- maybe select samples to contrast? color?
             deriving (Typeable,Data)

test, disp, info  :: Options
test = Test &= details ["Prints a test string.",
                       "This is useful for checking that your terminal supports the graphical characters needed to generate sparklines.  It should look like gradual transitions from one color to the next."]
disp = Disp &= details ["Show per sample consensus sequences.","Reads mpile-formatted input, and shows the consensus for each sample, one line each.  This can potentially result in very long lines, use 'samtools mpileup' with the '-r' option to restrict output."]
info = Info &= details ["Show expected information values.","This shows the information value for observing each allele, i.e. how diagnostic each site is between the two samples."]

main :: IO ()
main = do
  opts <- cmdArgsRun $ cmdArgsMode $ modes [test,disp,info]
          &= summary "Visualize information from 'samtools mpileup' as sparklines"
          &= program "sparks"
  inp <- BL.getContents
  let ms = map readPile1 $ BL.lines inp
  mapM_ putStrLn $ case opts of
     Test -> teststr
     Disp -> sparklines ms
     Info -> infoline ms

teststr :: [String]
teststr = [ concat [count2char [5,b,0,0] | b <- [0..10]]
          ++concat [count2char [0,5,b,0] | b <- [0..10]]
          ++ concat [count2char [0,0,5,b] |  b <- [0..10]]
          ++ count2char [0,0,0,10] ++ setSGRCode []]

infoline :: [MPileRecord] -> [String]
infoline ms = (cons:esivs)
  where
    cons = makeConsensus (IUPAC,1,5) ms
    esivs = map (concat . (++ rst)) $ transpose $ map esivstr ms
    rst = [setSGRCode []]

esivstr :: MPileRecord -> [String]
esivstr m = let
  (s1:s2:_) = counts m
  [t1,t2] = map covC [s1,s2]
  [maxpos1,maxpos2] = map (snd . last . sort) [zip xs [0..3] | xs <- [toList s1,toList s2]]
  es = [esiv1 1.64 0.01 (x,t1-x) (y,t2-y) | (x,y) <- zip (toList s1) (toList s2)]
  sorted = sort (zip es [0..3::Int])
  echar bg fg f = setSGRCode [SetColor Foreground Dull fg,SetColor Background Dull bg]
                  ++ if f<0 then [chr (0x2589-max 1 (min 8 (round (20*abs f))))]
                       else [chr (0x2580+max 1 (min 8 (round (10*f))))]
  toCol = ([Red,Blue,Green,Yellow]!!)
  in if sum (map abs es) > 0.1
     then let
       (plus,ppos) = last sorted
       (minus,mpos) = head sorted
       in [echar Black (toCol ppos) plus, echar (toCol mpos) Black minus]
     else [echar Black (toCol maxpos1) 0, echar (toCol maxpos2) Black 99]

sparklines :: [MPileRecord] -> [String]
sparklines = map sparkline . transpose . map counts

transpose :: [[a]] -> [[a]]
transpose xs
  | any null xs = []
  | otherwise = map head xs : transpose (map tail xs)

sparkline :: [Counts] -> String
sparkline = (++setSGRCode []) . concatMap (count2char . toList)

-- display a frequency spectrum as a sparkline char

count2char :: [Int] -> String
count2char [0,0,0,0] = setSGRCode [SetColor Background Dull Black] ++ " "
count2char cts = let
  ((a,b):(c,d):_) = reverse $ sort $ zip cts [0..3]
  [(majpos,major),(minpos,minor)] = sort [(b,a),(d,c)]
  freq = fromIntegral major / fromIntegral (major + minor)
  in f2char ("acgt"!!majpos) ("acgt"!!minpos) freq
     -- show (majpos,minpos,majfreq,counts) -- 

-- From colors and frequency, calculate the appropriate sparkline char
-- f should range be in the range [0..9]
f2char :: Char -> Char -> Double -> String
f2char c1 c2 f
  | f<=0 = setSGRCode [SetColor Foreground Dull Black,SetColor Background Dull bg] ++ [c2]
  | f>=1 = setSGRCode [SetColor Foreground Dull Black,SetColor Background Dull fg] ++ [c1]
  | otherwise = cols ++ [chr (0x2580+round (9*f))]
  where
    cols = setSGRCode [SetColor Foreground Dull fg,SetColor Background Dull bg]
    fg = toCol c1
    bg = toCol c2
    toCol 'a' = Red -- A
    toCol 'c' = Blue -- C
    toCol 'g' = Green -- G 
    toCol 't' = Yellow -- T
    toCol x = error ("illegal character "++show x)

reset :: IO ()
reset = setSGR []
