package algorithm

import (
	DNAType "GoDNA/algorithm/dnatype"
	"math"
	"math/rand"
)

type TSA struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
}

func NewTSA(maxIt int) Algorithm {
	return &TSA{
		Pop:          nil,
		MaxIteration: maxIt,
	}
}

func (tsa *TSA) GetName() string {
	return "TSA"
}

// Initialize implements Algorithm.
func (tsa *TSA) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	tsa.Pop = pop
}

// Iteration implements Algorithm.
func (tsa *TSA) Iteration(hyperPlaneNorm bool, cd bool) *DNAType.DNAAgent {
	fits := tsa.Pop.Fit()
	ZMin := tsa.Pop.ZMin()
	bestIndex, _ := NDKPSort(fits, ZMin, tsa.Pop.Size(), hyperPlaneNorm, cd)

	bestIndividual := tsa.Pop.At(bestIndex).Clone()
	gbest := bestIndividual.Variance()

	for it := range tsa.MaxIteration {
		oldPop := tsa.Pop.Clone()
		// xmin := 0.1
		// xmax := 1.5
		//gamma := 0.5 - 0.4*(1-(float64(it)/float64(tsa.MaxIteration)))
		gamma := 0.1 + 0.4*(1-float64(it)/float64(tsa.MaxIteration))
		popMean := mean(tsa.Pop)
		for i := range oldPop.Size() {
			x := oldPop.At(i).Clone().Variance()
			for j := range tsa.Pop.VarianceDim() {
				// A1 := (rand.Float64() + rand.Float64())
				// A1 -= rand.Float64() * rand.Float64()
				// xr := xmin + rand.Float64()*(xmax-xmin)
				// A1 /= xr
				A1 := C0(gamma)
				//c2 := rand.Float64()
				if rand.Float64() < 0.4 {
					c3 := rand.Float64()
					if c3 >= 0.5 {
						d_pos := math.Abs(x[j] - popMean[j]) //* c2
						x[j] = x[j] + A1*d_pos
					} else {
						d_pos := math.Abs(x[j] - popMean[j]) //* c2
						x[j] = x[j] - A1*d_pos
					}
				} else {
					c3 := rand.Float64()
					d_pos := math.Abs(x[j] - gbest[j])
					if c3 >= 0.5 {
						//d_pos := math.Abs(x[j] - c2*x[j])

						x[j] = x[j] + A1*d_pos
					} else {
						//d_pos := math.Abs(x[j] - c2*x[j])
						x[j] = x[j] - A1*d_pos
					}

				}
				// chosen := rand.Int() % oldPop.Size()
				// x[j] = x[j] + 0.5*(rand.Float64()-0.5)*tsa.Pop.At(chosen).Variance()[j]
				// x[j] = (x[j] + oldPop.At(chosen).Variance()[j]) / (1.5 + rand.Float64())

				x[j] = max(x[j], tsa.Pop.LB())
				x[j] = min(x[j], tsa.Pop.UB())
			}
			tsa.Pop.UpdatePosition(i, x)
		}
		tsa.Pop.PostWork()

		tsa.Pop.Join(oldPop)

		fits = tsa.Pop.Fit()
		ZMin = tsa.Pop.ZMin()
		var selectedIndex []int
		bestIndex, selectedIndex = NDKPSort(fits, ZMin, tsa.Pop.Size()/2, hyperPlaneNorm, cd)
		bestIndividual = tsa.Pop.At(bestIndex).Clone()
		gbest = bestIndividual.Variance()
		tsa.Pop.Select(selectedIndex)
	}
	return bestIndividual
}
