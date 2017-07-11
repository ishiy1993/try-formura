import Control.Monad
import Data.List
import System.Environment

-- gauss
x0 = 50
a = 10

-- sin
b0 = 1.0
b1 = 0.5
k = 2*pi/50

main :: IO ()
main = do
    [scheme] <- getArgs
    es <- forM [100*nx | nx <- [1..10]] $ \nx -> do
        let dx = 100 / fromIntegral nx
            xs = [(fromIntegral i)*dx | i <- [0..nx-1]]
            dds = derivatives scheme xs
            dds' = differences scheme dx xs
            error = dds |-| dds'
            l1 = foldl' (\acc (_,e) -> acc + abs e) 0 error
            nx' = fromIntegral nx
        return (dx,l1,l1/nx')
    writeData scheme es

type State = [(Double, Double)]

writeData :: String -> [(Double,Double,Double)] -> IO ()
writeData scheme es = do 
    let file = "data/" ++ scheme ++ "-derivative.err"
        st = unlines $ map (\(dx,l1,mre) -> show dx ++ " " ++ show l1 ++ " " ++ show mre) es
    writeFile file st

derivatives :: String -> [Double] -> State
derivatives scheme = map (derivative scheme)
    where
        derivative "gauss" x = (x, -2*(x-x0)*exp(-((x-x0)/a)**2) / a**2)
        derivative "sin" x = (x, -b1*k*cos(k*x))

differences :: String -> Double -> [Double] -> State
differences scheme dx xs = zipWith5 (\d__ d_ (x,_) d' d'' -> (x,(d__ - 8*d_ + 8*d' - d'')/(12*dx))) ds__ ds_ ds ds' ds''
    where
        ds = map (build scheme) xs
        ds__ = map snd $ rot (-2) ds
        ds_  = map snd $ rot (-1) ds
        ds'  = map snd $ rot 1 ds
        ds'' = map snd $ rot 2 ds

build "gauss" x = (x, 1 + exp(-((x-x0)/a)**2))
build "sin" x = (x, b0 + b1*sin(-k*x))

rot :: Int -> State -> State
rot i ds = as ++ bs
    where
        n = length ds
        j = if i >= 0 then i else n+i
        (bs,as) = splitAt j ds

(|-|) :: State -> State -> State
ds |-| ds' = zipWith (\(x,d) (x',d') -> (x,d-d')) ds ds'
