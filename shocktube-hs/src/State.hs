module State where

import qualified Data.Vector as V

data Cell = Cell
    { id :: Int
    , dens :: Double
    , mass :: Double
    , engy :: Double
    }

instance Show Cell where
    show (Cell id d m e) = unwords [show id, show d, show m, show e]

type State = V.Vector Cell

formatState :: State -> String
formatState = unlines . V.toList . V.map show
