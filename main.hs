{-  author: samtenka
 -  change: 2021-06-17
 -  create: 2021-06-17
 -  descrp: combinator-based parsing for the COW language 
 -  thanks: 
 -  to use: `make atom` 
 -}

import System.IO
import Data.List
import Data.Bits
import Data.Hashable
import qualified Data.HashSet as HSet
import Numeric.Natural

-- ============================================================================
-- ===  0. PERMUTATION GROUPS  ================================================
-- ============================================================================

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~~~~  0.0. Endomorphisms  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-- -------------------  0.0.0. an Endo is a endomap in {0...n-1}  -------------

newtype Endo = E (Natural, Natural -> Natural)

instance Show Endo where
    show (E(n,ff)) =
        foldr (\x-> \ss -> (show x)++ss) "" $ map ff [0..n-1] 

instance Eq Endo where
    (==) (E(n,ff)) (E(m,gg)) =
        (n==m) && (null differ)
        where differ = [k | k<-[0..n-1], (ff k) /= (gg k)] 

instance Ord Endo where
    (<=) (E(n,ff)) (E(m,gg)) =
        (n<=m) || ((n==m) && (null differ || (ff (head differ) <= gg (head differ))))
        where differ = [k | k<-[0..n-1], (ff k) /= (gg k)]

instance Hashable Endo where
    -- no need for security.  just need speed and improbable collisions
    hashWithSalt a (E(n,ff)) =
        sum [fromIntegral $ 101^k * (ff k) | k<-[0..n-1]]

card :: Endo -> Natural
card (E(n,ff)) = n

-- -------------------  0.0.1. we may apply and compose Endos  ----------------

appl :: Endo -> Natural -> Natural
ident :: Natural -> Endo
identFor :: Endo -> Endo
mult :: Endo -> Endo -> Endo
power :: Endo -> Natural -> Endo

appl (E(n,ff)) k = ff k

ident n = E(n, \a->a)

identFor (E(n,ff)) = ident n 

mult (E(n,ff)) (E(m,gg)) = E(n,ff.gg)

power f k = if k == 0 then identFor f
                      else if k `mod` 2 == 0 then double
                                             else mult f double  
                           where double = power (mult f f) (k `div` 2)

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~~~~  0.1. Permutations  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

isPerm :: Endo -> Bool
order :: Endo -> Natural
inverse :: Endo -> Endo
elOrder :: Endo -> Natural -> Natural
transposition :: Natural -> Natural -> Natural -> Endo

isPerm (E(n,ff)) = sort (map ff [0..n-1]) == [0..n-1]

order f = head $ [k | k <- [1..], power f k == identFor f]

inverse f = power f $ (order f) - 1

elOrder f x = head $ [k | k <- [1..], appl (power f k) x == x]

transposition n a b = E(n, \x -> if x==a then b
                                         else if x==b then a
                                                      else x)

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~~~~  0.2. Endomorphism Submonoids  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-- -------------------  0.2.0. an SMon is a set of endomorphisms  -------------

newtype SMon = S (HSet.HashSet Endo)

instance Show SMon where
    show (S(fs)) = (show $ HSet.size fs) ++ " : " ++ 
                   (HSet.foldr (\f-> \ss-> (show f)++" "++ss) "" fs)

size :: SMon -> Natural
size (S(fs)) = fromIntegral $ HSet.size fs

-- -------------------  0.2.1. submonoid constructors  ------------------------

trivial :: Natural -> SMon 
cyclic :: Endo -> SMon 
generate :: Natural -> [Endo] -> SMon 
lmul :: Endo -> SMon -> SMon 
rmul :: SMon -> Endo -> SMon 

trivial n = S(HSet.singleton $ ident n)

cyclic f = generate (card f) [f]

generate n fs =  
    S(grow fs (HSet.fromList fs) (HSet.singleton $ ident n))
    where grow fs frontier accum =
            if HSet.null frontier then accum
                                  else grow fs new_frontier new_accum
            where translates = [HSet.map (mult f) frontier | f<-fs]
                  new_accum = HSet.union frontier accum 
                  new_frontier = HSet.difference (foldr HSet.union (HSet.empty) translates) new_accum
            
lmul f (S(gs)) = S(HSet.map (mult f) gs) 

rmul (S(gs)) f = S(HSet.map (\g -> mult g f) gs) 

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~~~~  0.3. Permutation Groups  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

symmetric :: Natural -> SMon

-- we generate from a (connected, few-edged) *low-diameter* graph!
symmetric n = generate n [transposition n a b | a<-[0..0], b<-[1..n-1]] 

-- ============================================================================
-- ===  4. MAIN LOOP  =========================================================
-- ============================================================================

ex_a, ex_b :: Endo
ex_a = E(4, \k -> if k<4 then [1,0,3,4,2]!!(fromIntegral k) else k)
ex_b = E(4, \k -> if k<4 then [1,2,3,4,0]!!(fromIntegral k) else k)

sm_a = cyclic ex_a 
sm_b = cyclic ex_b 
sm_c = generate 5 [ex_a, ex_b] 
sm_d = symmetric 7
sm_e = symmetric 8

main = do putStr "hello!\n\n"
          putStr $ show sm_a ++ "\n\n"
          putStr $ show sm_b ++ "\n\n"
          putStr $ show sm_c ++ "\n\n"
          putStr $ (show $ size sm_d) ++ " elements \n\n"
          putStr $ (show $ size sm_e) ++ " elements \n\n"
          contents <- readFile "valediction"
          putStr contents
