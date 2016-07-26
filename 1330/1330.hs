import Control.Applicative
import Control.Monad
import Data.Array
import Data.Maybe

import qualified Data.ByteString.Char8 as B

int = fst . fromJust . B.readInt

readInt :: IO Int
readInt = int <$> B.getLine

readIntPair :: IO [Int]
readIntPair = parseLine <$> B.getLine
    where
        parseLine = (map int) . B.words

getAnswer prefixSums (from:to:_) = (prefixSums ! to) - (prefixSums ! (from - 1))

genPrefixSums numbers = listArray (0, length numbers) (scanl (+) 0 numbers)

main = do
    n <- readInt
    numbers <- replicateM n readInt
    let prefixSums = genPrefixSums numbers
    m <- readInt
    queries <- replicateM m readIntPair
    mapM_ (print . getAnswer prefixSums) queries

