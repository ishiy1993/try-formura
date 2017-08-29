import Data.List
import System.Environment

main :: IO ()
main = do
    [fn] <- getArgs
    let ps = getParams fn
    ds <- readData fn
    let res = map (calcViscosity ps) ds
    let (dir,base) = drop 1 <$> break (=='/') fn
        fn' = dir ++ "/viscosity-" ++ base
    writeData fn' res

data Params = Params
    { alpha :: Double
    , alpha' :: Double
    , dx :: Double
    }

data Info = Info
    { pos :: Double
    , idens :: Double
    , velc :: Double
    , pres :: Double
    , idensX :: Double
    , velcX :: Double
    , presX :: Double
    }

getParams :: FilePath -> Params
getParams fn = Params a a' h
    where
        split = unfoldr (\s -> if null s then Nothing else Just (drop 1 <$> break (=='-') s))
        ps = split fn
        a = read $ ps !! 3
        a' = read $ ps !! 4
        n = read $ ps !! 5
        h = 100/n

readData :: FilePath -> IO [Info]
readData fn = do
    ds <- map (map read . words) . lines <$> readFile fn
    return $ map toInfo ds

toInfo :: [Double] -> Info
toInfo l = Info
    { pos = l !! 0
    , idens = l !! 1
    , velc = l !! 3
    , pres = l !! 4
    , idensX = l !! 5
    , velcX = l !! 6
    , presX = l !! 7
    }

calcViscosity :: Params -> Info -> (Double, Double)
calcViscosity (Params a a' h) (Info x b u p bx ux px)
    | ux < 0 = (x, -a*c*h*ux/b + a'*h*h*ux*ux/b)
    | otherwise = (x, 0)
    where
        gm = 1.4
        c = (gm*b*p)**(1/2)

writeData :: FilePath -> [(Double,Double)] -> IO ()
writeData fn ds = writeFile fn $ unlines $ map (\(x,v) -> show x ++ " " ++ show v) ds
