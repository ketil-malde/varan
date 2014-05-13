-- High performance counting data structure

module Count (Counts(..)
             , addA, addC, addG, addT, addV
             , addA_, addC_, addG_, addT_                                       
             , getA, getC, getG, getT
             , covC, ptAdd, ptSum
             , toList, sumList 
             ) where

import Variants       
import Data.Int
import Data.List (foldl1')

data Counts = C { getA_, getC_, getG_, getT_ :: {-# UNPACK #-} !Int32
                , getV :: ![Variant]
                }

getA, getC, getG, getT :: Integral i => Counts -> i
getA = fromIntegral . getA_
getC = fromIntegral . getC_
getG = fromIntegral . getG_
getT = fromIntegral . getT_
{-# INLINE getA #-}
{-# INLINE getC #-}
{-# INLINE getG #-}
{-# INLINE getT #-}

addA, addC, addG, addT :: Integral i => Counts -> i -> Counts
addA c i = c { getA_ = getA_ c + fromIntegral i }
addC c i = c { getC_ = getC_ c + fromIntegral i }
addG c i = c { getG_ = getG_ c + fromIntegral i }
addT c i = c { getT_ = getT_ c + fromIntegral i }
{-# INLINE addA #-}
{-# INLINE addC #-}
{-# INLINE addG #-}
{-# INLINE addT #-}

addA_, addC_, addG_, addT_ :: Counts -> Int32 -> Counts
addA_ c i = c { getA_ = getA_ c + i }
addC_ c i = c { getC_ = getC_ c + i }
addG_ c i = c { getG_ = getG_ c + i }
addT_ c i = c { getT_ = getT_ c + i }


addV :: Counts -> Variant -> Counts
addV c v = c { getV = v : getV c }
{-# INLINE addV #-}

covC :: Counts -> Int
covC c = fromIntegral (getA_ c + getC_ c + getG_ c + getT_ c)

ptAdd :: Counts -> Counts -> Counts
ptAdd a b = a `addA` (getA_ b) `addC` (getC_ b) `addG` (getG_ b) `addT` (getT_ b)

ptSum :: [Counts] -> Counts
ptSum = foldl1' ptAdd

-- Convert 'Counts' to a list of allele counts
toList :: Num a => Counts -> [a]
toList x = map fromIntegral [getA_ x,getC_ x,getG_ x,getT_ x]

-- | Pointwise summation of the input lists
sumList :: Num a => [[a]] -> [a]
sumList = foldr (zipWith (+)) [0,0,0,0]
{-# DEPRECATED sumList "use ptSum" #-}

