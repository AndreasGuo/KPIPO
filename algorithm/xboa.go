package algorithm

import (
	DNAType "GoDNA/algorithm/dnatype"
	"math"
	"math/rand"
)

type XBOA struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
	c            float64
	a            float64
	p            float64
}

func (xboa *XBOA) GetName() string {
	return "XBOA"
}

func (xboa *XBOA) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	xboa.Pop = pop
	xboa.c = 0.02
	xboa.a = 0.2
	xboa.p = 0.8
}

func (xboa *XBOA) Iteration(hyperPlaneNorm bool, origin bool, cd bool) *DNAType.DNAAgent {
	// logger := log.Default()
	// islog := false
	fits := xboa.Pop.Fit()
	ZMin := xboa.Pop.ZMin()
	bestIndex, _ := NDKPSort(fits, ZMin, xboa.Pop.Size(), hyperPlaneNorm, cd)

	bestIndividual := xboa.Pop.At(bestIndex)

	for it := 0; it < int(xboa.MaxIteration); it++ {
		oldPop := xboa.Pop.Clone()
		extraIndex := oldPop.Size()
		for i := 0; i < int(oldPop.Size()); i++ {
			if i == bestIndex {
				continue
			}
			// frangrace by hm
			f := xboa.c * (math.Pow(fits[i][2], xboa.a))
			x := oldPop.At(i).Variance()
			if rand.Float64() > xboa.p {
				metaIndex := i
				for metaIndex == i && xboa.Pop.Size() > 1 {
					metaIndex = rand.Intn(oldPop.Size())
				}
				mate := oldPop.At(metaIndex).Variance()
				crossoverPoint := rand.Intn(len(x))
				x1 := make([]float64, crossoverPoint, len(x))
				copy(x1, x[:crossoverPoint])
				x2 := make([]float64, len(x)-crossoverPoint, len(x))
				copy(x2, x[crossoverPoint:])
				mate1 := make([]float64, crossoverPoint, len(x))
				mate2 := make([]float64, len(x)-crossoverPoint, len(x))
				copy(mate1, mate[:crossoverPoint])
				copy(mate2, mate[crossoverPoint:])

				x1 = append(x1, mate2...)
				mate1 = append(mate1, x2...)
				x2 = mate1

				// boundary control
				for j := 0; j < xboa.Pop.VarianceDim(); j++ {
					x1[j] = max(x1[j], xboa.Pop.LB())
					x2[j] = max(x2[j], xboa.Pop.LB())
					x1[j] = min(x1[j], xboa.Pop.UB())
					x2[j] = min(x2[j], xboa.Pop.UB())
				}

				xboa.Pop.UpdatePosition(i, x1)
				xboa.Pop.Append([]*DNAType.DNAAgent{&DNAType.DNAAgent{DNARowVariance: x2}})
				xboa.Pop.UpdatePosition(extraIndex, x2)
				extraIndex = 1
			} else {
				r1 := rand.Float64()
				r2 := rand.Float64()
				j := rand.Intn(oldPop.Size())
				k := rand.Intn(oldPop.Size())
				gamma := 0.5 - 0.4*(1-(float64(it)/float64(xboa.MaxIteration)))
				for d := range len(x) {
					x[d] += f * (r1*r2*oldPop.At(j).Variance()[d] - oldPop.At(k).Variance()[d])
					if rand.Float64() < 0.1 {
						x[d] += C0(gamma)
					}
				}
				for j := 0; j < xboa.Pop.VarianceDim(); j++ {
					x[j] = max(x[j], xboa.Pop.LB())
					x[j] = min(x[j], xboa.Pop.UB())
				}
				xboa.Pop.UpdatePosition(i, x)
			}
		}
		xboa.c += 0.025 / (xboa.c + float64(xboa.MaxIteration))
		xboa.Pop.PostWork()
		xboa.Pop.Join(oldPop)

		fits = xboa.Pop.Fit()
		ZMin = xboa.Pop.ZMin()
		bestIndex, selectedIndex := NDKPSort(fits, ZMin, oldPop.Size(), hyperPlaneNorm, cd)
		bestIndividual = xboa.Pop.At(bestIndex)
		xboa.Pop.Select(selectedIndex)
	}
	return bestIndividual
}
