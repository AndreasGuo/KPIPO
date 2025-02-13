package algorithm

import (
	DNAType "GoDNA/algorithm/dnatype"
	"log"
	"math"
	"math/rand"
)

type PO struct {
	Pop          *DNAType.DNAPopulation
	MaxIteration int
}

func NewPO(maxIteration int) Algorithm {
	po := &PO{Pop: nil, MaxIteration: maxIteration}
	return po
}

func (po *PO) GetName() string {
	return "PO"
}

func (po *PO) Initialize(pop *DNAType.DNAPopulation, inds ...*DNAType.DNAAgent) {
	pop.Init()
	if len(inds) > 0 {
		pop.Append(inds)
	}
	po.Pop = pop
}

// PO + NDSort + Knee Point
func (po *PO) Iteration(hyperPlaneNorm bool, cd bool) *DNAType.DNAAgent {
	logger := log.Default()
	islog := false
	fits := po.Pop.Fit()
	ZMin := po.Pop.ZMin()
	bestIndex, _ := NDKPSort(fits, ZMin, po.Pop.Size(), hyperPlaneNorm, cd)
	origin := false

	bestIndividual := po.Pop.At(bestIndex).Clone()
	gbest := po.Pop.At(bestIndex).Variance()
	sita := rand.Float64() * math.Pi

	for it := 0; it < int(po.MaxIteration); it++ {
		// fmt.Print("\r", it)
		oldPop := po.Pop.Clone()
		popMean := mean(po.Pop)
		for i := 0; i < int(po.Pop.Size()); i++ {
			// if i == selectedIndex {
			// 	continue
			// }
			st := rand.Intn(4)
			x := po.Pop.At(i).Variance()
			if islog {
				logger.Println("it: ", it, "st: ", st)
				logger.Println("origin variance: ", x)
			}

			switch st {
			case 0:
				if origin {
					ost0(x, gbest, popMean, po.Pop.VarianceDim(), it, po.MaxIteration)
				} else {
					st0(x, gbest, popMean, po.Pop.VarianceDim(), it, po.MaxIteration)
				}

			case 1:
				if origin {
					ost1(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration)
				} else {
					st1(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration)
				}
			case 2:
				if origin {
					ost2(x, popMean, po.Pop.VarianceDim(), it, po.MaxIteration)
				} else {
					st2(x, popMean, po.Pop.VarianceDim(), it, po.MaxIteration)
				}
			case 3:
				if origin {
					ost3(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration, sita)
				} else {
					st3(x, gbest, po.Pop.VarianceDim(), it, po.MaxIteration, sita)
				}
			}
			out := make([]int, len(x))
			for k := range x {
				out[k] = int(math.Round(x[k]))
			}
			if islog {
				logger.Println("After st: ", out)
			}
			//boundary control
			for j := 0; j < po.Pop.VarianceDim(); j++ {
				x[j] = max(x[j], po.Pop.LB())
				x[j] = min(x[j], po.Pop.UB())
			}
			for k := range x {
				out[k] = int(math.Round(x[k]))
			}
			if islog {
				logger.Println("After boundary: ", out)
			}
			po.Pop.UpdatePosition(i, x)
		}
		po.Pop.PostWork()

		po.Pop.Join(oldPop)

		fits = po.Pop.Fit()
		ZMin = po.Pop.ZMin()
		var selectedIndex []int
		bestIndex, selectedIndex = NDKPSort(fits, ZMin, po.Pop.Size()/2, hyperPlaneNorm, cd)
		bestIndividual = po.Pop.At(bestIndex).Clone()
		gbest = po.Pop.At(bestIndex).Variance()
		po.Pop.Select(selectedIndex)

		sita = rand.Float64() * math.Pi
	}

	return bestIndividual //po.Pop.At(selectedIndex)
}

func levy(lh int) []float64 {
	beta := 1.5
	sigma := math.Gamma(1+beta) * math.Sin(math.Pi*beta/2)
	sigma /= math.Gamma(1+beta/2) * beta * math.Pow(1, (beta-1)/2)
	sigma = math.Pow(sigma, 1/beta)
	u := make([]float64, lh)
	for i := range u {
		u[i] = sigma * rand.NormFloat64() / math.Pow(math.Abs(rand.NormFloat64()), 1/beta)
	}

	//for i := range u {
	//	u[i] /= 1.1
	//	//u[i] *= 1.2
	//}
	return u
}

func ost0(x, gbest, popMean []float64, dim, it, maxIt int) {
	levyDim := levy(dim)
	for i := 0; i < dim; i++ {
		r := rand.Float64()
		x[i] = (x[i] - gbest[i]) * levyDim[i] * r
		x[i] += (1 - r) * popMean[i] * math.Pow(1-float64(it)/float64(maxIt), 2*float64(it)/float64(maxIt))
	}
}

func st0(x, gbest, popMean []float64, dim, it, maxIt int) {
	levyDim := levy(dim)
	p1 := float64(4)
	for i := 0; i < dim; i++ {
		r := rand.Float64()
		x[i] += (x[i] - gbest[i]) * levyDim[i] * r / p1
		x[i] += (1 - r) * popMean[i] * math.Pow(1-float64(it)/float64(maxIt), 2*float64(it)/float64(maxIt)) / p1
	}
}

func ost1(x, gbest []float64, dim, it, maxIt int) {
	levyDim := levy(dim)
	for i := 0; i < dim; i++ {
		r := rand.NormFloat64()
		x[i] += gbest[i]*levyDim[i] + r
	}
}

func st1(x, gbest []float64, dim, it, maxIt int) {
	levyDim := levy(dim)
	p2 := float64(2)
	for i := 0; i < dim; i++ {
		r := rand.NormFloat64()
		x[i] += gbest[i]*levyDim[i]/p2 + r*(1-float64(it)/float64(maxIt))
	}
}

func ost2(x, popMean []float64, dim, it, maxIt int) {
	p := rand.Float64()
	alpha := rand.NormFloat64() / 5
	if p <= 0.5 {
		for i := 0; i < dim; i++ {
			x[i] += alpha * (1 - float64(it)/float64(maxIt)) * (x[i] - popMean[i])
		}
	} else {
		for i := 0; i < dim; i++ {
			x[i] += alpha * math.Exp(float64(-it)/(rand.Float64()*float64(maxIt)))
		}
	}
}

func st2(x, popMean []float64, dim, it, maxIt int) {
	p := rand.Float64()
	alpha := rand.NormFloat64() / 2
	if p <= 0.5 {
		for i := 0; i < dim; i++ {
			x[i] += alpha * (1 - float64(it)/float64(maxIt)) * (x[i] - popMean[i])
		}
	} else {
		for i := 0; i < dim; i++ {
			x[i] += alpha * math.Exp(float64(-it)/(rand.Float64()*float64(maxIt))) * (x[i] - popMean[i])
		}
	}
}

func ost3(x, gbest []float64, dim, it, maxIt int, sita float64) {
	pow := math.Pow(float64(it)/float64(maxIt), 2/float64(maxIt))
	for i := 0; i < dim; i++ {
		r := rand.Float64()
		x[i] += r*math.Cos(math.Pi*float64(it)/(2*float64(maxIt)))*(gbest[i]-x[i]) -
			math.Cos(sita)*(x[i]-gbest[i])*pow
	}
}

func st3(x, gbest []float64, dim, it, maxIt int, sita float64) {
	pow := math.Pow(1-(float64(it)/float64(maxIt)), 2/float64(maxIt))
	for i := 0; i < dim; i++ {
		r := rand.Float64()
		x[i] += r*math.Cos(math.Pi*float64(it)/(2*float64(maxIt)))*(gbest[i]-x[i]) -
			math.Cos(sita)*(x[i]-gbest[i])*pow
	}
}

func mean(pop *DNAType.DNAPopulation) []float64 {
	size := pop.Size()
	dim := pop.VarianceDim()
	out := make([]float64, dim)
	for i := range size {
		v := pop.At(i).Variance()
		for j := range dim {
			out[j] += v[j]
		}
	}

	for j := range dim {
		out[j] /= float64(size)
	}

	return out
}
