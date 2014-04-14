-- High performance counting data structure

module Count (Counts(..), getcounts
             , getA, getC, getG, getT
             , addA, addC, addG, addT
             , covC, ptAdd, ptSum
             , toList, sumList 
             , Int64) where

import Variants       
       
import Data.Int
import Data.Bits
import Data.List (foldl1')

data Counts = C {-# UNPACK #-} !Int64 [Variant]

getcounts :: Counts -> Int64
getcounts (C cs _) = cs

getA, getC, getG, getT :: Int64 -> Int
getA = flip shiftR 48 .fromIntegral . (.&.) 0xFFFF000000000000 
getC = flip shiftR 32 . fromIntegral . (.&.) 0xFFFF00000000 
getG = flip shiftR 16 . fromIntegral . (.&.) 0xFFFF0000 
getT = fromIntegral . (.&.) 0xFFFF

addA, addC, addG, addT :: Int64 -> Int -> Int64
addA c i = c + fromIntegral (i `shiftL` 48)
addC c i = c + fromIntegral (i `shiftL` 32)
addG c i = c + fromIntegral (i `shiftL` 16)
addT c i = c + fromIntegral i

covC :: Int64 -> Int
covC c = getA c + getC c + getG c + getT c

ptAdd :: Int64 -> Int64 -> Int64
ptAdd a b = a `addA` (getA b) `addC` (getC b) `addG` (getG b) `addT` (getT b)

ptSum :: [Int64] -> Int64
ptSum = foldl1' ptAdd

-- Convert 'Counts' to a list of allele counts
toList :: Num a => Int64 -> [a]
toList x = map fromIntegral [getA x,getC x,getG x,getT x]
{-# DEPRECATED toList "use covC and sumC" #-}

-- | Pointwise summation of the input lists
sumList :: Num a => [[a]] -> [a]
sumList = foldr (zipWith (+)) [0,0,0,0]
{-# DEPRECATED sumList "use covC instead" #-}

