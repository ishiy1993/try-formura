{-# LANGUAGE DuplicateRecordFields #-}
module Main where

import Control.Monad (foldM_, when)
import qualified Data.Vector as V
import Text.Printf

import Config
import State

initialState :: V.Vector Cell
initialState = V.generate (xR-xL+1) $ \i ->
    if i <= (xR+xL) `div` 2
       then Cell i densL massL engyL
       else Cell i densR massR engyR

main :: IO ()
main = do
    let ss = iterate step (0, initialState)
    save 300 ss $ \i s -> do
        when (i `mod` 10 == 0) $ do
            let t = fst s
            putStrLn $ printf "i = %d, t = %f" i t
            let f = printf "data/%f.dat" t
            writeFile f $ formatState s
        return $ i+1

type Vec = (Double, Double, Double)

(|+|) :: Vec -> Vec -> Vec
(a1,a2,a3) |+| (b1,b2,b3) = (a1+b1,a2+b2,a3+b3)
infixl 6 |+|

(|-|) :: Vec -> Vec -> Vec
(a1,a2,a3) |-| (b1,b2,b3) = (a1-b1,a2-b2,a3-b3)
infixl 6 |-|

(*|) :: Double -> Vec -> Vec
c *| (b1,b2,b3) = (c*b1,c*b2,c*b3)
infixl 7 *| 

(|/) :: Vec -> Double -> Vec
(a1,a2,a3) |/ c = (a1/c,a2/c,a3/c)
infixl 7 |/ 

toVec :: Cell -> Vec
toVec (Cell _ d m e) = (d,m,e)

flux :: Vec -> Vec
flux (d,m,e) = (dF,mF,eF)
    where
        dF = m
        mF = (gamma - 1)*e + (3 - gamma)*m**2/(2*d)
        eF = gamma*e*m/d - (gamma - 1)*m**3/(2*d**2)

maxV :: V.Vector Cell -> Double
maxV = V.maximum . V.map calcV
    where
        calcV :: Cell -> Double
        calcV (Cell _ d m e) = let u = m/d
                                   p = (gamma-1)*(e - m**2/(2*d))
                                in abs u + sqrt (gamma * p / d)

step :: State -> State
step = lax

lax :: State -> State
lax (t0, ss) = (t1, V.map laxScheme ss)
    where
        dt = 0.4 * dx / maxV ss
        t1 = t0 + dt
        laxScheme (Cell i _ _ _)
            | i == xL = Cell i densL massL engyL
            | i == xR = Cell i densR massR engyR
            | otherwise =
                let qL = toVec $ ss V.! (i-1)
                    qR = toVec $ ss V.! (i+1)
                    fL = flux qL
                    fR = flux qR
                    (d,m,e) = (qL |+| qR)|/2 |-| (dt/dx)*|(fR |-| fL)|/2
                 in Cell i d m e

fds :: State -> State
fds = undefined

save :: Int -> [State] -> (Int -> State -> IO Int) -> IO ()
save n ss act = foldM_ act 0 $ take n ss
