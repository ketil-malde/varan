module AgrestiCoull where

-- zscores: 1.654 is one-sided 95%
--          1.96 is two-sided
--          2.326 is one-sided 99%

confidenceInterval :: Double -> Int -> Int -> (Double,Double)
confidenceInterval z succs fails = let
  s = fromIntegral succs
  f = fromIntegral fails
  z2 = z * z
  nest = s + f + z2
  pest = (s + z2/2)/nest
  delta = z*sqrt(pest*(1-pest)/nest)
  in (pest-delta,pest+delta)