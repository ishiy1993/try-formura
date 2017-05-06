{-# LANGUAGE DuplicateRecordFields #-}
module Main where

import Control.Monad (foldM_, when)
import Data.Maybe
import qualified Data.Vector as V
import System.Environment (getArgs)
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
    as <- getArgs
    let scheme = fromMaybe "lax" $ listToMaybe as
    putStrLn scheme
    let ss = iterate (step scheme) (0, initialState)
    save 300 ss $ \i s -> do
        when (i `mod` 10 == 0) $ do
            let t = fst s
            putStrLn $ printf "i = %d, t = %f" i t
            let f = printf "data/%s-%f.dat" scheme t
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

dot :: Vec -> Vec -> Double
dot (a1,a2,a3) (b1,b2,b3) = a1*b1 + a2*b2 + a3*b3

toVec :: Cell -> Vec
toVec (Cell _ d m e) = (d,m,e)

flux :: Vec -> Vec
flux (d,m,e) = (dF,mF,eF)
    where
        dF = m
        mF = (gamma - 1)*e + (3 - gamma)*m**2/(2*d)
        eF = gamma*e*m/d - (gamma - 1)*m**3/(2*d**2)

type Matrix = V.Vector Double

mkMatrix :: Vec -> Vec -> Matrix
mkMatrix (dL,mL,eL) (dR,mR,eR) = r >*< l >*< r'
    where
        l = V.fromList [abs (u_-c_), 0, 0
                       , 0, abs u_, 0
                       , 0, 0, abs (u_+c_)]
        r = V.fromList [ 1, 1, 1
                       , u_-c_, u_, u_+c_
                       , h_-u_*c_, u_**2/2, h_+u_*c_]
        r' = V.fromList [ (b1+u_/c_)/2, -(1/c_+b2*u_)/2, b2/2
                        , 1-b1, b2*u_, -b2
                        , (b1-u_/c_)/2, (1/c_-b2*u_)/2, b2/2]
        b1 = u_**2/2 * (gamma-1)/c_**2
        b2 = (gamma-1)/c_**2
        dL' = sqrt dL
        dR' = sqrt dR
        d_ = dL' * dR'
        u_ = (uL*dL' + uR*dR')/(dL' + dR')
        c_ = sqrt $ (gamma-1)*(h_ - u_**2/2)
        h_ = (hL*dL' + hR*dR')/(dL' + dR')
        uL = mL/dL
        uR = mR/dR
        pL = (gamma-1)*(eL-mL**2/2/dL)
        pR = (gamma-1)*(eR-mR**2/2/dR)
        hL = (eL+pL)/dL
        hR = (eR+pR)/dR

(>!) :: Matrix -> (Int,Int) -> Double
m >! (i,j) = m V.! (3*i+j)

(>*<) :: Matrix -> Matrix -> Matrix
m1 >*< m2 = V.generate 9 $ \n ->
    let i = n `div` 3
        j = n `mod` 3
     in sum $ map (\k -> m1 >! (i,k) * m2 >! (k,j)) [0..2]

(>*|) :: Matrix -> Vec -> Vec
m >*| (a0,a1,a2) = (b0, b1, b2)
    where b0 = m >! (0,0) * a0 + m >! (0,1) * a1 + m >! (0,2) * a2
          b1 = m >! (1,0) * a0 + m >! (1,1) * a1 + m >! (1,2) * a2
          b2 = m >! (2,0) * a0 + m >! (2,1) * a1 + m >! (2,2) * a2

maxV :: V.Vector Cell -> Double
maxV = V.maximum . V.map calcV
    where
        calcV :: Cell -> Double
        calcV (Cell _ d m e) = let u = m/d
                                   p = (gamma-1)*(e - m**2/(2*d))
                                in abs u + sqrt (gamma * p / d)

step :: String -> State -> State
step "lax" = lax
step "fds" = fds
step "muscl" = muscl

lax :: State -> State
lax (t0, ss) = (t1, V.map next ss)
    where
        dt = 0.4 * dx / maxV ss
        t1 = t0 + dt
        next (Cell i _ _ _)
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
fds (t0, ss) = (t1, V.map next ss)
    where
        dt = 0.4 * dx / maxV ss
        t1 = t0 + dt
        next (Cell i d m e)
            | i == xL = Cell i densL massL engyL
            | i == xR = Cell i densR massR engyR
            | otherwise =
                let q = (d,m,e)
                    qL = toVec $ ss V.! (i-1)
                    qR = toVec $ ss V.! (i+1)
                    f = flux q
                    fL = flux qL
                    fR = flux qR
                    aL = mkMatrix qL q
                    aR = mkMatrix q qR
                    fL' = (f |+| fL |-| aL >*| (q |-| qL))|/2
                    fR' = (fR |+| f |-| aR >*| (qR |-| q))|/2
                    (d',m',e') = q |-| (dt/dx) *| (fR' |-| fL')
                 in Cell i d' m' e'

muscl :: State -> State
muscl (t0,ss) = (t1, next ss)
    where
        dt = 0.4 * dx / maxV ss
        t1 = t0 + dt
        sub = V.zipWith (\(Cell i d m e) (fd,fm,fe) -> Cell i (d-fd) (m-fm) (e-fe))
        next st = let st1 = sub st $ df (dt/dx/4) st
                      st2 = sub st $ df (dt/dx/3) st1
                      st3 = sub st $ df (dt/dx/2) st2
                  in  sub st $ df (dt/dx) st3
            where
                df a s = V.map (\(Cell i d m e) -> a *| dfs i) s
                    where
                        dfs i = f i |-| f (i-1)
                        q i | i <= xL = (densL,massL,engyL)
                            | i >= xR = (densR,massR,engyR)
                            | otherwise = toVec $ s V.! i
                        f i = let dq = q i |-| q (i-1)
                                  dQ = q (i+1) |-| q i
                                  ep = 10**(-6)
                                  s = (2*dot dq dQ + ep)/(dot dq dq + dot dQ dQ + ep)
                                  qL = q i |+| ((s/4) *| ((1-s/3)*|dq |+| (1+s/3)*|dQ))
                                  qR = q (i+1) |-| ((s/4) *| ((1-s/3)*|dQ |+| (1+s/3)*|dq))
                                  a = mkMatrix qL qR
                              in (flux qR |+| flux qL |-| a >*| (qR |-| qL))|/2

save :: Int -> [State] -> (Int -> State -> IO Int) -> IO ()
save n ss act = foldM_ act 0 $ take n ss
