import qualified Data.Vector as V
import System.Environment (getArgs)
import Text.Printf

import Config

velcL = massL/densL
velcR = massR/densR
presL = (gamma-1)*(engyL - densL*velcL**2/2)
presR = (gamma-1)*(engyR - densR*velcR**2/2)
cL = sqrt (gamma*presL/densL)
cR = sqrt (gamma*presR/densR)

f :: Double -> Double
f p = sqrt (2/(gamma*(gamma-1+(gamma+1)*p))) * (p-1)
    - 2/(gamma-1) * cL/cR * (1-(presR*p/presL)**((gamma-1)/(2*gamma)))
    - (velcR-velcL)/cR

df :: Double -> Double
df p = sqrt (2/gamma*(gamma-1+(gamma+1)*p))
     * (1 - (gamma+1)*(p-1)/(2*(gamma-1+(gamma+1)*p)))
     + 2/(gamma-1) * a * (presR*p/presL)**a / p
    where a = (gamma-1)/(2*gamma)

solveP :: Double -> Double
solveP p | abs (p'-p) / p < 10**(-6) = p
         | otherwise = solveP p'
    where p' = p - f p / df p

main :: IO ()
main = do
    as <- getArgs
    case as of
      [] -> putStrLn "Need given time"
      (t:_) -> do
          let es = exactSolution (read t)
              fn = "data/exact-" ++ t ++ ".dat"
          writeFile fn $ formatState es

exactSolution :: Double -> State
exactSolution t = V.generate (xR-xL+1) $ build t p
    where p = solveP 3

build :: Double -> Double -> Int -> Basic
build t p i
    | i' < x0 + (velcL-cL)*t = Basic i dens5 velc5 pres5
    | i' <= x0 + v'*t = Basic i dens4 velc4 pres4
    | i' <= x0 + vcd*t = Basic i dens3 velc3 pres3
    | i' <= x0 + vs*t = Basic i dens2 velc2 pres2
    | otherwise = Basic i dens1 velc1 pres1
    where
        i' = fromIntegral i
        x0 = fromIntegral (xR+xL)/2
        v' = (gamma+1)*vcd/2 - cL - (gamma-1)*velcL/2
        vcd = velc2
        vs = velcR + (p-1)*cR**2/(gamma*(velc2-velcR))
        c4 = cL - (gamma-1)*(velc4 - velcL)/2
        dens1 = densR
        velc1 = velcR
        pres1 = presR
        dens2 = densR*(gamma-1+(gamma+1)*p)/(gamma+1+(gamma-1)*p)
        velc2 = velcR + cR*(p-1)*sqrt (2/(gamma*(gamma-1+(gamma+1)*p)))
        pres2 = presR*p
        dens3 = densL*(pres3/presL)**(1/gamma)
        velc3 = velc2
        pres3 = pres2
        dens4 = densL*(pres4/presL)**(1/gamma)
        velc4 = 2/(gamma+1) * ((i'-x0)/t + cL + (gamma-1)*cL/2)
        pres4 = presL*(c4/cL)**(2*gamma/(gamma-1))
        dens5 = densL
        velc5 = velcL
        pres5 = presL


data Basic = Basic
    { id :: Int
    , dens :: Double
    , velc :: Double
    , pres :: Double
    }

instance Show Basic where
    show (Basic i d v p) = printf "%d %f %f %f" i d v p

type State = V.Vector Basic

formatState :: State -> String
formatState = unlines . V.toList . V.map show
