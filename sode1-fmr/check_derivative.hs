import Control.Monad
import Data.List
import Numeric.AD
import System.Environment

-- gauss
x0 :: Floating a => a
x0 = 50.0
a :: Floating a => a
a = 10.0

gauss :: Floating a => a -> a
gauss x = 1 + exp(-((x-x0)/a)**2)

gauss' :: Floating a => a -> a
gauss' x = -2*(x-x0)*exp(-((x-x0)/a)**2) / a**2

gauss'' :: Floating a => a -> a
gauss'' = diff (diff gauss)

gauss''' :: Floating a => a -> a
gauss''' = diff (diff gauss')
-- gauss''' = diff (diff (diff gauss)) -- This result is NaN

-- sin
b0 :: Floating a => a
b0 = 1.0
b1 :: Floating a => a
b1 = 0.5
k :: Floating a => a
k = 2*pi/50

sinWave :: Floating a => a -> a
sinWave x = b0 + b1*sin(-k*x)

sinWave' :: Floating a => a -> a
sinWave' x = -b1*k*cos(k*x)

sinWave'' :: Floating a => a -> a
sinWave'' = diff (diff sinWave)

sinWave''' :: Floating a => a -> a
sinWave''' = diff (diff (diff sinWave))


fun :: Floating a => String -> a -> a
fun "gauss" = gauss
fun _ = sinWave

fun' :: Floating a => String -> a -> a
fun' "gauss" = gauss'
fun' _ = sinWave'

fun'' :: Floating a => String -> a -> a
fun'' "gauss" = gauss''
fun'' _ = sinWave''

fun''' :: Floating a => String -> a -> a
fun''' "gauss" = gauss'''
fun''' _ = sinWave'''


main :: IO ()
main = do
    [scheme] <- getArgs
    let f = fun scheme
    let f' = fun' scheme
    let f'' = fun'' scheme
    let f''' = fun''' scheme
    es <- forM [100*nx | nx <- [1..10]] $ \nx -> do
        let dx = 100 / fromIntegral nx
            xs = [(fromIntegral i)*dx | i <- [0..nx-1]]
            ys = map f xs
            ys' = map f' xs
            d2s = map f'' xs
            d2s' = diff2 dx ys ys'
            d3s = map f''' xs
            d3s' = diff3 dx ys'
            es2 = zipWith (-) d2s d2s'
            es3 = zipWith (-) d3s d3s'
            d2L1 = foldl' (\acc e -> acc + abs e) 0 es2
            d3L1 = foldl' (\acc e -> acc + abs e) 0 es3
            nx' = fromIntegral nx
        return (dx,d2L1/nx',d3L1/nx')
    writeData scheme es

writeData :: (Show a, Floating a) => String -> [(a,a,a)] -> IO ()
writeData scheme es = do 
    let file = "data/" ++ scheme ++ "-derivative.err"
        st = unlines $ map (\(dx,d2,d3) -> show dx ++ " " ++ show d2 ++ " " ++ show d3) es
    writeFile file st

diff2 :: Floating a => a -> [a] -> [a] -> [a]
diff2 dx ys ys' = zipWith5 (\yB y yA y'B y'A -> 2*(yA + yB - 2*y - dx*(y'A - y'B)/4)/dx/dx) ysB ys ysA ys'B ys'A
    where
        ysA = tail ys ++ [head ys]
        ysB = last ys : init ys
        ys'A = tail ys' ++ [head ys']
        ys'B = last ys' : init ys'

diff3 :: Floating a => a -> [a] -> [a]
diff3 dx ys' = zipWith3 (\yB y yA -> (yA + yB - 2*y)/dx/dx) ys'B ys' ys'A
    where
        ys'A = tail ys' ++ [head ys']
        ys'B = last ys' : init ys'
