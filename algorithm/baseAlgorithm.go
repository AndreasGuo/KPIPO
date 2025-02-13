package algorithm

import DNAType "GoDNA/algorithm/dnatype"

type anyFuncType func(...interface{}) any

type Algorithm interface {
	// Initialize a population, with fitness calculated
	// The fitness function, it takes all individual as input and calculate fitness parallel
	// and returns the size*objs matrix of fitness values.
	Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent)
	Iteration(hyperPlaneNorm bool, cd bool) *DNAType.DNAAgent
	GetName() string
}

//type Population[T int | float64] struct {
//	Size        int
//	Lb, Ub      T
//	VarianceDim int
//	ObjectDim   int
//	Pop         []Individual
//}

// type Population interface {
// 	Size() int
// 	Init()
// 	Append([]Individual)
// 	Fit() [][]float64
// 	ZMin() []float64
// 	At(int) Individual
// 	VarianceDim() int
// 	ObjectiveDim() int
// 	UpdatePosition(int, []float64)
// 	LB() float64
// 	UB() float64
// 	PostWork()
// 	Variance() [][]float64
// 	Clone() Population
// 	Join(Population)
// 	Select([]int)
// }

// type Individual interface {
// 	// Objs Returns object list
// 	Objs() []float64
// 	Variance() []float64
// 	UpdatePosition([]float64)
// 	// repair and fixGC, hook function
// 	PostWork()
// 	// this for dna Seq
// 	Represent() []int
// 	// this for dna seq tostr()
// 	String() (string, error)
// 	Clone() Individual
// 	SetObjs([]float64)
// }
