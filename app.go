package main

import (
	"GoDNA/algorithm"
	DNAType "GoDNA/algorithm/dnatype"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"gonum.org/v1/gonum/stat"
)

var zmin []float64

var startTime time.Time
var logFileName string

func App(config Config) {
	startTime = time.Now()
	// var alg = algorithm.NewTSA(config.MAXIT)
	var alg = algorithm.NewPO(config.MAXIT)
	logFileName = startTime.Format("2006-01-02=15=04") + alg.GetName() + ".txt"
	file, err := createFile(logFileName)
	if err != nil {
		fmt.Println("error on creating file: " + err.Error())
		return
	}
	defer file.Close()

	// create working pool
	fitChan := DNAType.CreateWorker(100, 100, 10)
	defer fitChan.Close()

	// init dna set
	var dnaSet = randomDNASet(config.DIM, config.DNASIZE)

	//global zmin, which is the min value of each objective dimension
	zmin = make([]float64, 5)
	for j := range zmin {
		zmin[j] = 400
	}

	result := printDNASet(dnaSet, fitChan)
	file.Write([]byte(result))

	for it := 0; it < config.DNASETITERATION; it++ {
		fmt.Println("DNASet iteration ", it+1, "/", config.DNASETITERATION)
		file.Write([]byte("DNASet iteration " + strconv.Itoa(it+1) + "/" + strconv.Itoa(config.DNASETITERATION) + ", "))
		// init objs
		// cause one sequence is changed
		// the others must be recauculated
		for i := range dnaSet {
			fitFunc := DNAType.FitnessCall(dnaSet, i, fitChan, config.MINVALUE, config.FITREVERSE)
			singleInvSlice := []*DNAType.DNAAgent{dnaSet[i]}
			fits, _ := fitFunc(singleInvSlice)
			dnaSet[i].SetObjs(fits[0])
		}

		// choose a sequence in DNA set to optimize
		index := chooseInvToOpt(dnaSet, config.MINVALUE)
		fmt.Println("To Optimize: ", index)
		file.Write([]byte("To Optimize: " + strconv.Itoa(index) + "\n\n"))

		fitFunc := DNAType.FitnessCall(dnaSet, index, fitChan, config.MINVALUE, config.FITREVERSE)
		// alg := algorithm.NewPO(config.MAXIT)
		//var alg algorithm.Algorithm = &algorithm.XBOA{Pop: nil, MaxIteration: config.MAXIT}
		//var alg = algorithm.NewGWO(config.MAXIT)
		pop := new(DNAType.DNAPopulation)
		pop.SetConfig(config.POPSIZE-1, config.DIM, 5, float64(config.LB), float64(config.UB))
		pop.SetFitFunc(fitFunc)
		alg.Initialize(pop, dnaSet[index])
		inv := alg.Iteration(config.PLANENORM, config.CD)
		dnaSet[index] = inv

		// generation log
		result := printDNASet(dnaSet, fitChan)
		file.Write([]byte(result))
		if it == config.DNASETITERATION-1 {
			endTime := time.Now()
			runningDuration := endTime.Sub(startTime).Seconds()
			result := "fit_reversed=" + strconv.FormatBool(config.FITREVERSE) + "\n"
			result += "planenorm=" + strconv.FormatBool(config.PLANENORM) + "\n"
			result += "DNA_set_iteration=" + strconv.Itoa(config.DNASETITERATION) + "\n"
			result += "pop_iteration=" + strconv.Itoa(config.MAXIT) + "\n"
			result += "pop_size=" + strconv.Itoa(config.POPSIZE) + "\n"
			result += "crowding_dis=" + strconv.FormatBool(config.CD) + "\n"
			// result += "original PO =" + strconv.FormatBool(config.ORIGINPO) + "\n"
			result += "running time (s)=" + strconv.FormatFloat(runningDuration, 'f', 4, 64)
			file.Write([]byte(result))
			//saveResult(result, alg.GetName())
		}
	}
}

func createFile(filename string) (*os.File, error) {
	resultsDir := "results"
	if _, ok := os.Stat(resultsDir); os.IsNotExist(ok) {
		err := os.Mkdir(resultsDir, os.ModePerm)
		if err != nil {
			fmt.Println("error on creating resutls dir: ", err)
			return nil, err
		}
	}
	return os.Create(filepath.Join(resultsDir, filename))
}

// func saveResult(result string, alg string) {
// 	//now := time.Now()
// 	resultsDir := "results"
// 	if _, ok := os.Stat(resultsDir); os.IsNotExist(ok) {
// 		err := os.Mkdir(resultsDir, os.ModePerm)
// 		if err != nil {
// 			fmt.Println("error on creating resutls dir: ", err)
// 			return
// 		}
// 	}
// 	fileName := filepath.Join(resultsDir, str+".txt")
// 	err := os.WriteFile(fileName, []byte(result), 0644)
// 	if err != nil {
// 		fmt.Println("error on writing result: ", err)
// 		fmt.Println("result is:", result)
// 	}
// }

func printDNASet(dnaSet []*DNAType.DNAAgent, fitChan *DNAType.FitChan) string {
	var result = ""
	var setCt, setHp, setHm, setSm float64 = 0, 0, 0, 0
	var mts = make([]float64, 0, 7)
	var setMT float64 = 0
	for ind, inv := range dnaSet {
		DNAString, err := inv.String()
		if err != nil {
			panic("error while decoding")
		}
		result += strconv.Itoa(ind) + " " + DNAString

		// continuity
		fitChan.CtIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		ct := (<-fitChan.CtRe).Value
		result += " " + strconv.FormatFloat(ct, 'f', 4, 64)
		setCt += ct

		// hairpin
		fitChan.HpIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		hp := (<-fitChan.HpRe).Value
		result += " " + strconv.FormatFloat(hp, 'f', 4, 64)
		setHp += hp

		// hm
		hmList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			fitChan.HmIn <- DNAType.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
			hmList[j] = (<-fitChan.HmRe).Value
		}
		hm := sum(hmList)
		result += " " + strconv.FormatFloat(hm, 'f', 4, 64)
		setHm += hm

		// sm
		smList := make([]float64, len(dnaSet))
		for j, o := range dnaSet {
			if j != ind {
				fitChan.SmIn <- DNAType.SeqMapPair{Index1: ind, Index2: j, Seq1: inv.Represent(), Seq2: o.Represent()}
				smList[j] = (<-fitChan.SmRe).Value
			}
		}
		sm := sum(smList)
		result += " " + strconv.FormatFloat(sm, 'f', 4, 64)
		setSm += sm

		//mt
		fitChan.MtIn <- DNAType.SeqMapSingle{Index: ind, Seq: inv.Represent()}
		mt := (<-fitChan.MtRe).Value
		mts = append(mts, mt)
		setMT += mt
		result += " " + strconv.FormatFloat(mt, 'f', 4, 64) + "\n"
	}

	for range 2 + len((*dnaSet[0]).Seq) {
		result += " "
	}
	result += " " + strconv.FormatFloat(setCt/7.0, 'f', 4, 64)
	result += " " + strconv.FormatFloat(setHp/7.0, 'f', 4, 64)
	result += " " + strconv.FormatFloat(setHm/7.0, 'f', 4, 64)
	result += " " + strconv.FormatFloat(setSm/7.0, 'f', 4, 64)
	result += " " + strconv.FormatFloat(setMT/7.0, 'f', 4, 64)
	result += "(" + strconv.FormatFloat(stat.Variance(mts, nil), 'f', 4, 64) + ")"
	result += "\n"
	fmt.Println(result)
	return result
}

func randomDNASet(dim, size int) []*DNAType.DNAAgent {
	dnaSet := []*DNAType.DNAAgent{}
	for i := 0; i < size; i++ {
		dnaSet = append(dnaSet, DNAType.CreateDNAAgent(dim, 0, 3))
	}
	return dnaSet
}

func sum[T int | float64](lt []T) T {
	var s T = 0
	for i := range lt {
		s += lt[i]
	}
	return s
}

func chooseInvToOpt(dnaSet []*DNAType.DNAAgent, minVal float64) (idx int) {
	if len(dnaSet) == 0 {
		panic("no inv")
	}

	//zmax := make([]float64, len(zmin))

	for _, dna := range dnaSet {
		objs := dna.Objs()
		for i := range len(zmin) {
			zmin[i] = min(zmin[i], objs[i])
			//zmax[i] = max(zmax[i], objs[i])
		}
	}
	fmt.Println("zmin: ", zmin)

	distance := make([]float64, len(dnaSet))
	for i, dna := range dnaSet {
		objs := dna.Objs()
		distance[i] = 1
		for j := range zmin {
			//distance[i] += math.Pow(objs[j]-zmin[j], 2)
			// if zmax[j] == zmin[j] {
			// 	continue
			// }
			distance[i] *= max(1, (objs[j] - zmin[j]))
		}
	}
	fmt.Println("distance: ", distance)

	dix := 0.0
	for i := range distance {
		if distance[i] > dix {
			idx = i
			dix = distance[i]
		}
	}
	return idx
}

// func variance(nums []float64) float64 {
// 	var sum float64 = 0
// 	for _, n := range nums {
// 		sum += n
// 	}
// 	mean := sum / float64(len(nums))
// 	var v float64 = 0
// 	for _, n := range nums {
// 		v += (n - mean) * (n - mean)
// 	}
// 	return v / (float64(len(nums) - 1))
// }
