package main

import (
	"flag"
)

func main() {
	// flags
	var dnaLength = flag.Int("vdim", 20, "length of sigle DNA sequence")
	var dnaSetSize = flag.Int("setsize", 7, "size of DNA set")
	var popSize = flag.Int("popsize", 50, "size of population")
	var maxIt = flag.Int("maxit", 200, "max iteration each algorithm run")
	var dnaSetIteration = flag.Int("optit", 30, "number of optimize iterations")
	var minVal = flag.Float64("minval", 2e-2, "min value when calculate product items which need avoid zero")
	var fitReverse = flag.Bool("fitrev", true, "is reversed src and dst in calculating hm and sm")
	var planeNorm = flag.Bool("norm", false, "normalize objs before calculate distance between individuals and hyperplane")
	var chooseToOpt = flag.Int("worstdef", 0, "the method to choose worst sequence in DNA set, 0-product; 1-elucid distance")
	// var originPO = flag.Bool("originpo", false, "if use origin po")
	var cd = flag.Bool("cd", false, "whether use crowding distance instead of knee point")
	flag.Parse()

	// boundary of problem
	// DNA bases are encoded into int number
	// which are 0,1,2,3
	const lb = 0
	const ub = 3

	config := Config{
		DIM:             *dnaLength,
		DNASIZE:         *dnaSetSize,
		POPSIZE:         *popSize,
		MAXIT:           *maxIt,
		LB:              lb,
		UB:              ub,
		DNASETITERATION: *dnaSetIteration,
		MINVALUE:        *minVal,
		FITREVERSE:      *fitReverse,
		PLANENORM:       *planeNorm,
		CHOOSETOOPT:     *chooseToOpt,
		// ORIGINPO:        *originPO,
		CD: *cd,
	}

	App(config)

}
