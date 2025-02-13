package algorithm

import (
	DNAType "GoDNA/algorithm/dnatype"
	"math"
	"math/rand"
)

type OTSA struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
}

func NewOTSA(maxIt int) Algorithm {
	return &OTSA{
		Pop:          nil,
		MaxIteration: maxIt,
	}
}

func (tsa *OTSA) GetName() string {
	return "OTSA"
}

// Initialize implements Algorithm.
func (tsa *OTSA) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	tsa.Pop = pop
}

// Iteration implements Algorithm.
func (tsa *OTSA) Iteration(hyperPlaneNorm bool, cd bool) *DNAType.DNAAgent {
	fits := tsa.Pop.Fit()
	ZMin := tsa.Pop.ZMin()
	bestIndex, _ := NDKPSort(fits, ZMin, tsa.Pop.Size(), hyperPlaneNorm, cd)

	bestIndividual := tsa.Pop.At(bestIndex).Clone()
	gbest := bestIndividual.Variance()

	for range tsa.MaxIteration {
		oldPop := tsa.Pop.Clone()
		xmin := 1.0
		xmax := 4.0
		for i := range oldPop.Size() {
			x := oldPop.At(i).Clone().Variance()
			for j := range tsa.Pop.VarianceDim() {
				A1 := (rand.Float64() + rand.Float64())
				A1 -= rand.Float64() * rand.Float64()
				xr := xmin + rand.Float64()*(xmax-xmin)
				A1 /= xr
				c2 := rand.Float64()
				c3 := rand.Float64()

				d_pos := math.Abs(x[j] - gbest[j]*c2)
				if c3 >= 0.5 {
					x[j] = x[j] + A1*d_pos
				} else {
					x[j] = x[j] - A1*d_pos
				}

				chosen := rand.Int() % oldPop.Size()
				x[j] = (x[j] + oldPop.At(chosen).Variance()[j]) / (2 + rand.Float64())

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
