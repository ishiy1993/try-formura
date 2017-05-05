module State where

import qualified Data.Vector as V
import Text.Printf

import Config (gamma)

data Cell = Cell
    { id :: Int
    , dens :: Double
    , mass :: Double
    , engy :: Double
    }

instance Show Cell where
    show (Cell id d m e) = printf "%d %f %f %f %f %f %f" id d m e v p c
        where v = m/d
              p = (gamma-1)*(e - d*v**2/2)
              c = sqrt $ gamma*p/d

type State = (Double, V.Vector Cell)

formatState :: State -> String
formatState = unlines . V.toList . V.map show . snd
