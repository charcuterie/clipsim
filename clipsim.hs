import Control.Monad
import Data.Array.ST
import Data.Array.Unboxed
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as M
import Data.Maybe
import Numeric.Statistics (pvar)
import Options.Applicative
import System.Directory
import System.FilePath

----------------------------------------
-- Command-line argument parsing code --
----------------------------------------

data Options = Options
  { chrom    :: String   -- interval chromosome
  , cStart   :: Int      -- interval start site
  , cEnd     :: Int      -- interval end site
  , kernel   :: Bool     -- whether to apply RBF kernel
  , fileList :: String } -- list of input files

opts :: ParserInfo Options
opts = info (helper <*> options)
  ( fullDesc
 <> progDesc ("Calculate the pairwise similarity of clip-peak data over " ++
              "the interval CHROM:M-N and print to STDOUT" )
 <> header "clipsim - a tool for comparing CLIP-peak data" )

options :: Parser Options
options = Options
  <$> chromOption
  <*> startOption
  <*> endOption
  <*> kernelOption
  <*> inputOption

inputOption :: Parser String
inputOption = strOption
              ( long "input"
             <> short 'i'
             <> metavar "FILE"
             <> help "List of files to compare; one filename per line" )

chromOption :: Parser String
chromOption = strOption
              ( long "chrom"
             <> metavar "CHROM"
             <> help ("The chromosome of the region of interest, " ++ 
                     "e.g., \"chr2\"") )

startOption :: Parser Int
startOption = option auto
              ( long "start"
             <> metavar "M"
             <> help "The start coordinate of the region of interest" )

endOption :: Parser Int
endOption = option auto
            ( long "end"
           <> metavar "N"
           <> help "The end coordinate of the region of interest" )

kernelOption :: Parser Bool
kernelOption = switch
               ( long "kernel"
              <> help "Apply the radial-basis function kernel to the output" )

-------------------------
-- BED-processing code --
-------------------------

-- BED intervals are a chrom, start coord, and end coord
type BedInterval = (String, Int, Int)

-- Parse a BED file to a list of BED intervals
parseBed :: String -> [BedInterval]
parseBed = map (getCoords . words) . lines
  where
    getCoords (a:b:c:_) = (a, read b, read c)
    getCoords _         = error "Malformed BED file" -- die immediately (?)

-- Intersection of a BED interval with a list of BED intervals
intersectBeds :: BedInterval -> [BedInterval] -> [BedInterval]
intersectBeds b = catMaybes . map (intersectBed b)

-- Get the intersection of two BED intervals
intersectBed :: BedInterval -> BedInterval -> Maybe BedInterval
intersectBed (chr0, st0, end0) (chr1, st1, end1)
  | chr0 /= chr1               = Nothing
  | st0 >= end1 || st1 >= end0 = Nothing
  | otherwise                  = Just $ (chr0, (max st0 st1), (min end0 end1))

-- Do these two BED intervals overlap?
overlaps :: BedInterval -> BedInterval -> Bool
overlaps x y = intersectBed x y /= Nothing

-- Convert a list of BED intervals into an uarray of booleans, true if the
-- position is within at least one interval. The first argument is the
-- interval of interest, used to set the dimension of the uarray.
getCoverage :: BedInterval -> [BedInterval] -> UArray Int Bool
getCoverage i@(iChr, iSt, iEnd) bs = runSTUArray $ do
  arr <- newArray (iSt, iEnd) False
  forM_ bs $ \b@(bChr, bSt, bEnd) ->
    when (b `overlaps` i) $
      forM_ [(max iSt bSt)..(min iEnd bEnd)] $ \k ->
        writeArray arr k True
  return arr

---------------------------
-- Pairwise-scoring code --
---------------------------

type ScoreMethod = Bool -> Bool -> Double

hamming :: ScoreMethod
hamming i j
  | i /= j    = 1
  | otherwise = 0


score :: ScoreMethod                    -- the scoring method
      -> Map FilePath (UArray Int Bool) -- a mapping from filepath to cov array
      -> (FilePath, FilePath)           -- a pairing of filepaths
      -> Double                         -- the score
score s m (f0, f1) = sum $ (zipWith s) covs0 covs1
  where
    covs0 = entries $ M.lookup f0 m
    covs1 = entries $ M.lookup f1 m
    entries Nothing  = []
    entries (Just x) = elems x

-----------------
-- Kernel code --
-----------------

data Kernel = RbfKernel | NoKernel

applyKernel :: Kernel -> [Double] -> [Double]
applyKernel NoKernel  xs = xs
applyKernel RbfKernel xs = map (applyKernel' var) xs
  where
    var = pvar xs
    applyKernel' var x = exp 1 ** (negate $ x ** 2 / (2 * var))

--------------------------
-- Printing/output code --
--------------------------

printScores :: [(FilePath, FilePath)] -> [Double] -> IO ()
printScores ((f0, f1):fs) (x:xs) = do
  putStrLn $ (takeBaseName f0) ++ "\t" ++ (takeBaseName f1) ++ "\t" ++ (show x)
  printScores fs xs
printScores _ _ = return ()

-------------------
-- Program entry --
-------------------

{-  Line-by-line summary
 -
 -  1.  Parse the command-line arguments.
 -  2.  Read the filepaths from the input text file.
 -  3.  Make the "interval-of-interest" (ioi), i.e., the gene to calculate over
 -  4.  For each filepath, make a list of bed records that overlap the ioi.
 -  5.  For each bed-record list, make an boolean array, indicating peaks
 -  6.  Map each filepath to its peak array
 -  7.  Generate every pairwise combination of filepaths
 -  8.  Select a kernel.
 -  9.  Score each filepath-pair and apply the kernel.
 -  10. Print the results.
 -}

main :: IO ()
main = do
  args <- execParser opts
  paths <- (readFile . fileList) args >>= return . lines
  let ioi = (chrom args, cStart args, cEnd args)
  beds <- mapM (liftM ((intersectBeds ioi) . parseBed) . readFile) paths
  let covs = map (getCoverage ioi) beds
  let path2cov = M.fromList $ zip paths covs
  let pathPairs = [(p0, p1) | p0 <- paths, p1 <- paths, p0 < p1]
  let kern = if (kernel args) then RbfKernel else NoKernel
  let scores = applyKernel kern $ map (score hamming path2cov) pathPairs
  printScores pathPairs scores
