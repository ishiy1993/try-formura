module Main where

import Data.List (intersperse)
import Data.Vector (Vector, generate, (!))
import qualified Data.Vector as V

main :: IO ()
main = do
  let xs0 = buildInitialState 30
      res = take 10 $ iterate (step 0.5) xs0
  writeFile "out.dat" $ format res

type State = Vector (Int, Double)

buildInitialState :: Int -> State
buildInitialState x0 = generate 100 (\n -> if n < x0 then (n,1) else (n,0))

step :: Double -> State -> State
step nu xs = V.map (\(i,u) -> (i, u - nu*(uR i - uL i)/2)) xs
    where
        uR 99 = 0
        uR i = snd $ xs ! (i+1)
        uL 0 = 1
        uL i = snd $ xs ! (i-1)

format :: [State] -> String
format = unlines . intersperse "" . map (unlines . map (\(i,u) -> show i ++ " " ++ show u) . V.toList)
