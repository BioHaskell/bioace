{- |
   Read ACE format assembly files

   These are typically output by sequence assembly tools,
   like CAP3 or Phrap.

   Documented in the section labelled \"ACE FILE FORMAT\" at
   <http://bozeman.mbt.washington.edu/consed/distributions/README.14.0.txt>

   Briefly: each field is a line starting with a two letter code,
            in some cases followed by data lines termintated by a blank line.
   Here's an brief example how an ACE file looks like:

   @
          AS contigs reads
          CO contig_name bases reads segments compl (CAP3: segments=0)
          sequence
          BQ base_qualities
          AF read1 compl padded_start_consensus (negatives meaning?)
          AF read2 ..
          BS segments
          RD read1 bases info_items info_tags (latter two set to 0 by CAP3)
          sequence
          QA read1 qual_start qual_end align_start align_end
          DS (phred header? left empty by CAP3)
          RD read2 ...
   @

   As far as I know, this is only used for nucleotide sequences.
-}

{-# LANGUAGE CPP #-}

module Bio.Alignment.Ace (readACE, writeACE, Assembly(..), ptest, reads ) where

import Prelude   hiding (lines,words,readFile,unwords -- ByteString clashes
                        ,reads,pred                   -- Assembly clash
                        )

import Bio.Core.Sequence
import Bio.Core.Strand
import Bio.Alignment.AlignData (Sequence(..), Gaps, Alignment, extractGaps)

import qualified Data.ByteString.Lazy.Char8 as B
import Data.ByteString.Lazy.Char8 (ByteString,words,pack,unpack,readFile,unwords)

import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Pos (newPos)

import Control.Monad (liftM) -- ,when?
import Data.Char (chr)

instance Show Sequence where
  show (Seq lab seq qual) = unpack $ unSD seq

data Assembly = Asm { contig :: (Sequence, Gaps), fragments :: Alignment}
                deriving Show

{-# DEPRECATED reads "Stupid name, replaced by 'fragments'." #-}
reads :: Assembly -> Alignment
reads = fragments -- deprecated, stupid, stupid name.

type Str = ByteString

-- | ACE header lines with parameters
--   The tokenizer (scanner) should convert input into a list of these,
--   which in turn can be parsed by Parsec
data ACE = AS Str Str | CO Str Str Str Str Str
         | BQ | AF Str Str Str | BS [Str]
         | RD Str Str Str Str
         | QA Str Str Str Str Str
         | DS [Str] | Other [Str] | Empty
           deriving (Eq)

instance Show ACE where
    show (AS x y)       = "AS "++uw [x,y]
    show (CO a b c d e) = "CO "++uw [a,b,c,d,e]
    show BQ = "BQ"
    show (AF a b c) = "AF "++uw [a,b,c]
    show (RD a b c d) = "RD "++uw [a,b,c,d]
    show (QA a b c d e) = "QA "++uw [a,b,c,d,e]
    show (DS ss) = "DS "++uw ss
    show (Other ss) = uw ss
    show Empty = "(blank)"
    show _ = "unknown ACE string"

uw :: [Str] -> String
uw = unpack . unwords

-- | The Parsec parser type
type AceParser a = GenParser (SourcePos,ACE) () a

-- | Parse a single token, primitive parser
parse1 :: (ACE -> Maybe a) -> AceParser a
parse1 p = token sho pos pred
    where sho  (_,t) = show t
          pos  (n,_) = n
          pred (_,t) = p t

-- | Test parser p on a list of ACE elements
ptest :: Show a => String -> AceParser a -> [ACE] -> IO ()
ptest m p = parseTest p . source m

-- | Add SourcePoses to a stream of ACEs.
source :: String -> [ACE] -> [(SourcePos,ACE)]
source m = zip (iterate (\sp -> incSourceLine sp 1) (newPos m 1 0))

-- | Parse a complete ACE file as a set of assemblies.
ace :: AceParser [[Assembly]]
ace = many1 ace1

ace1 :: AceParser [Assembly]
ace1 = do
  as
  many blank
  many (do ctg >>= asm) -- apparently, CAP3 outputs empty assemblies (AS 0 0)

-- | parse the initial header
as :: AceParser (Int,Int)
as = parse1 (\t -> case t of AS cs rs -> do c <- B.readInt cs
                                            r <- B.readInt rs
                                            return (fst c,fst r)
                             _ -> Nothing) <?> "AS <int> <int>"

blank :: AceParser ()
blank = parse1 (\t -> case t of Empty -> Just ()
                                _ -> Nothing) <?> "empty line"

-- | parse the contig and quality information (CO, BQ)
ctg :: AceParser (Sequence, Gaps)
ctg = do
  name <- co
  sd   <- sdata
  let (sd',gaps) = extractGaps $ SeqData {unSD = sd}
  many blank
  -- Vector NTI produces ACE without BQ for the contig
  msq <- do bq
            sq <- qdata
            return (Just $ QualData {unQD = sq})
       <|> return Nothing
  many blank
  -- todo: gaps?
  return (Seq (SeqLabel {unSL = name}) sd' msq, gaps)

co, sdata, qdata :: AceParser Str
co = parse1 (\t -> case t of CO name a b c _compl -> do
                               _bs  <- B.readInt a
                               _rds <- B.readInt b
                               _seg <- B.readInt c
                               return name
                             _ -> Nothing) <?> "CO name <int> <int> <int> bool"

sdata = do return . B.concat =<< many1 sdata1
  where sdata1 = parse1 (\t -> case t of Other sd -> Just (unwords sd)
                                         _        -> Nothing) <?> "sequence data"

qdata = do return . B.concat =<< many1 qdata1
  where qdata1 = parse1 (\t -> case t of Other sd -> liftM (pack . map chr) (readInts sd)
                                         _        -> Nothing) <?> "sequence data"

-- | Read a list of Ints in the Maybe monad
readInts :: [ByteString] -> Maybe [Int]
readInts [] = Just []
readInts (x:xs) = do (i,_) <- B.readInt x
                     is    <- readInts xs
                     return $ (i:is)

bq :: AceParser ()
bq = parse1 (\t -> case t of BQ -> Just (); _ -> Nothing) <?> "BQ"

-- | Given the CO info, get the AFS'es
asm :: (Sequence, Gaps) -> AceParser Assembly
asm cg = do
  many blank
  afs <- many1 af
  _bss <- many bs
  many blank
  rds cg afs

-- | Parse a list of AFS, followed by actual read, and merge them
-- afs :: Sequence -> AceParser [Sequence] -- plus some auxiliary info?
af :: AceParser (Str, Strand, Offset)
af = parse1 (\t -> case t of AF a b c -> mkAF a b c
                             _        -> Nothing) <?> "AF name (U|C) pad_start"
    where mkAF a b c = do
            b' <- case unpack b of "U" -> Just Plus; "C" -> Just Minus; _ -> Nothing
            c' <- liftM (fromIntegral . fst) (B.readInt c)
            return (a,b',c')

bs :: AceParser (Int,Int,Str)
bs = parse1 (\t -> case t of BS [x,y,n] -> do 
                               x' <- readInt' x
                               y' <- readInt' y
                               return (x',y',n)
                             _     -> Nothing) <?> "BS x y name"

readInt' :: Str -> Maybe Int
readInt' = liftM (fromIntegral . fst) . B.readInt

rds :: (Sequence,Gaps) -> [(Str, Strand, Offset)] -> AceParser Assembly
rds cg xs = do
     r <- many1 rseq
     -- todo: check the number and merge with the afs
     let f (_name,d,off) (s,gs) = (off,d,s,gs)
     return $ Asm { contig = cg, fragments = zipWith f xs r }

rseq :: AceParser (Sequence, Gaps)
rseq = do
  (rn,_len,_,_) <- rd
  (s,gaps) <- return . (\x -> extractGaps $ SeqData {unSD = x}) =<< sdata
  -- when (B.length s == fromIntegral len) $ (fail "Incorrect sequence length!")
  -- todo: fix gaps!
  many1 blank
  qa
  ds
  many blank
  return (Seq (SeqLabel {unSL = rn}) s Nothing,gaps) -- huh?

-- | parse each read (RD, QA, DS)
--   Vector NTI appears to insert solitary RDs, sometimes even without any sequence data!?
--   This is not supported at this point.
rd :: AceParser (Str,Int,Int,Int)
rd = parse1 (\t -> case t of RD a b c d -> do [x,y,z] <- readInts [b,c,d]
                                              return (a,x,y,z)
                             _ -> Nothing) <?> "RD <string> <int> <int> <int>"

qa :: AceParser ()
qa = parse1 (\_ -> Just ())

ds :: AceParser ()
ds = parse1 (\t -> case t of DS _ -> Just (); _ -> Nothing) <?> "DS"

-- ----------------------------------------------------------
-- Convert lines into tokens
tokenize :: ByteString -> [ACE]
tokenize = map tokenize1 . B.lines

-- Tokenise a single line
-- todo: error on incorrect (partial) format, error reports with line number
tokenize1 :: ByteString -> ACE
tokenize1 l = case words l of
                [] -> Empty
                (h:ws) -> case (unpack h,ws) of
                            ("AS",[cs,rs])               -> AS cs rs
                            ("CO",[nm,bss,rs,segs,comp]) -> CO nm bss rs segs comp
                            ("BQ",[])                    -> BQ
                            ("AF",[a,b,c])             -> AF a b c
                            ("BS",_)                     -> BS ws
                            ("RD",[a,b,c,d])             -> RD a b c d
                            ("QA",[a,b,c,d,e])           -> QA a b c d e
                            ("DS",_)                     -> DS ws
                            _ -> Other (h:ws)

-- | Reading an ACE file.
readACE :: FilePath -> IO [[Assembly]]
readACE f = parseit =<< B.readFile f
    where parseit = \s -> case (parse ace f . source f . tokenize) s of
                            Left e  -> fail (show e)
                            Right a -> return a

-- formatError msg = error ("readACE: incorrect format in "++msg)

writeACE :: FilePath -> [Assembly] -> IO ()
writeACE = undefined

-- todo: hWrite etc

