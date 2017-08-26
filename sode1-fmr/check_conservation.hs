import Data.List (foldl', isPrefixOf, sortOn, unfoldr)
import System.Directory
import System.Environment

main :: IO ()
main = do
    [profile] <- getArgs
    let dir = "data"
    fs <- sortOn getTimeStep . map ((dir++"/")++) . filter (profile `isPrefixOf`) <$> listDirectory dir
    mapM_ process fs

getTimeStep :: FilePath -> Int
getTimeStep = read . fst . break (=='.') . last . split '-'
    where split c = unfoldr (\str -> if null str then Nothing else Just (drop 1 <$> break (==c) str))

process :: FilePath -> IO ()
process fn = do
    ds <- readData fn
    let (m,e) = foldl' calc (0,0) ds
    putStrLn $ unwords [fn, show m, show e]

readData :: FilePath -> IO [Info]
readData fn = do
    ds <- map (map read . words) . lines <$> readFile fn
    return $ map toInfo ds

calc :: (Double, Double) -> Info -> (Double, Double)
calc (m0, e0) d = (m0+m, e0+e)
    where
        gm = 1.4
        m = (*) <$> dens <*> velc $ d
        e = (pres d)/(gm-1) + (dens d)*(velc d)**2/2

data Info = Info
    { posX :: Double
    , dens :: Double
    , velc :: Double
    , pres :: Double
    }

toInfo :: [Double] -> Info
toInfo l = Info x d u p
    where
        x = l !! 0
        d = 1 / (l !! 1)
        u = l !! 2
        p = l !! 3
