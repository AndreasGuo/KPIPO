package DNAType

// author: Andreas Guo
// This file applied Repair function
// which make sure GC Content fixed at 50%
// Repair function requires length of DNA sequence is an even number

import "math/rand"

func count(seq Seq, condition func(int) bool) int {
	sum := 0
	for _, base := range seq {
		if condition(base) {
			sum += 1
		}
	}
	return sum
}

func isCorG(base int) bool {
	if base == C || base == G {
		return true
	}
	return false
}

func isAorT(base int) bool {
	if base == A || base == T {
		return true
	}
	return false
}

func Repair(seq Seq) Seq {
	countCG := count(seq, isCorG)
	countAT := count(seq, isAorT)

	if countCG > countAT {
		repairLen := (countCG - countAT) / 2
		for i := 0; i < repairLen; i++ {
			for {
				position := rand.Intn(len(seq))
				if seq[position] == C || seq[position] == G {
					if rand.Float32() > 0.5 {
						seq[position] = A
					} else {
						seq[position] = T
					}
					break
				}
			}

		}
	}

	if countCG < countAT {
		repairLen := (countAT - countCG) / 2
		for i := 0; i < repairLen; i++ {
			for {
				position := rand.Intn(len(seq))
				if seq[position] == A || seq[position] == T {
					if rand.Float32() > 0.5 {
						seq[position] = C
					} else {
						seq[position] = G
					}
					break
				}
			}
		}
	}

	return seq
}
