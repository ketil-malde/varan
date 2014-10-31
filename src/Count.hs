-- High performance counting data structure

module Count (Counts(..)
             , addA, addC, addG, addT, addN, addDel, addV
             , addA_, addC_, addG_, addT_, addN_, addDel_
             , getA, getC, getG, getT, getN, getDel
             , czero
             , covC, ptAdd, ptSum
             , toList, sumList 
             ) where

import Variants       
import Data.Int
import Data.List (foldl1')

data Counts = C { getA_, getC_, getG_, getT_, getN_, getDel_ :: {-# UNPACK #-} !Int32
                , getV :: ![Variant]
                }
czero :: Counts
czero = C 0 0 0 0 0 0 []

getA, getC, getG, getT, getN, getDel :: Integral i => Counts -> i
getA = fromIntegral . getA_
getC = fromIntegral . getC_
getG = fromIntegral . getG_
getT = fromIntegral . getT_
getN = fromIntegral . getN_
getDel = fromIntegral . getDel_
{-# INLINE getA #-}
{-# INLINE getC #-}
{-# INLINE getG #-}
{-# INLINE getT #-}
{-# INLINE getN #-}
{-# INLINE getDel #-}

addA, addC, addG, addT, addN, addDel :: Integral i => Counts -> i -> Counts
addA c i = c { getA_ = getA_ c + fromIntegral i }
addC c i = c { getC_ = getC_ c + fromIntegral i }
addG c i = c { getG_ = getG_ c + fromIntegral i }
addT c i = c { getT_ = getT_ c + fromIntegral i }
addN c i = c { getN_ = getN_ c + fromIntegral i }
addDel c i = c { getDel_ = getDel_ c + fromIntegral i }
{-# INLINE addA #-}
{-# INLINE addC #-}
{-# INLINE addG #-}
{-# INLINE addT #-}
{-# INLINE addN #-}
{-# INLINE addDel #-}

addA_, addC_, addG_, addT_, addN_, addDel_ :: Counts -> Int32 -> Counts
addA_ c i = c { getA_ = getA_ c + i }
addC_ c i = c { getC_ = getC_ c + i }
addG_ c i = c { getG_ = getG_ c + i }
addT_ c i = c { getT_ = getT_ c + i }
addN_ c i = c { getN_ = getN_ c + i }
addDel_ c i = c { getDel_ = getDel_ c + i }

addV :: Counts -> Variant -> Counts
addV c v = c { getV = v : getV c }
{-# INLINE addV #-}

covC :: Counts -> Int
covC c = fromIntegral (getA_ c + getC_ c + getG_ c + getT_ c) -- + getN_ c + getDel_ c)
 -- mostly, we want to compare the identified bases.

ptAdd :: Counts -> Counts -> Counts
ptAdd a b = a `addA` (getA_ b) `addC` (getC_ b) `addG` (getG_ b) `addT` (getT_ b) `addN` (getN_ b) `addDel` (getDel_ b) -- todo: concat variants?

ptSum :: [Counts] -> Counts
ptSum = foldl1' ptAdd

-- Convert 'Counts' to a list of allele counts
toList :: Num a => Counts -> [a]
toList x = map fromIntegral [getA_ x,getC_ x,getG_ x,getT_ x]

-- | Pointwise summation of the input lists
sumList :: Num a => [[a]] -> [a]
sumList = foldr (zipWith (+)) [0,0,0,0]
{-# DEPRECATED sumList "use ptSum" #-}

