package DNAType

import (
	"math"
)

type DNAPopulation struct {
	size         int
	varianceDim  int
	objectiveDim int
	lb, ub       float64
	individuals  []*DNAAgent
	fit          FitFuncType
	zmin         []float64
}

// select index in [indicies] of the pop
func (pop *DNAPopulation) Select(indicies []int) {
	pop.size = len(indicies)
	newIndividuals := make([]*DNAAgent, 0)
	for i := range indicies {
		newIndividuals = append(newIndividuals, pop.At(i))
	}
	pop.individuals = newIndividuals
}

func (pop *DNAPopulation) Clone() *DNAPopulation {
	newPop := DNAPopulation{
		size:         pop.size,
		varianceDim:  pop.varianceDim,
		objectiveDim: pop.objectiveDim,
		lb:           pop.lb,
		ub:           pop.ub,
		fit:          pop.fit,
		zmin:         pop.zmin,
	}
	newIndividuals := make([]*DNAAgent, pop.size)
	for i := range newIndividuals {
		newIndividuals[i] = pop.individuals[i].Clone()
	}
	newPop.individuals = newIndividuals
	return &newPop
}

func (pop *DNAPopulation) Join(population *DNAPopulation) {
	newIndividuals := make([]*DNAAgent, population.Size())
	for i := range population.Size() {
		newIndividuals[i] = population.At(i).Clone()
	}
	pop.Append(newIndividuals)
	pop.size = len(pop.individuals)
}

func (pop *DNAPopulation) Variance() [][]float64 {
	out := make([][]float64, pop.Size())
	for i := range len(pop.individuals) {
		v := make([]float64, pop.varianceDim)
		copy(v, pop.individuals[i].Variance())
		out[i] = v
	}
	return out
}

func (pop *DNAPopulation) SetFitFunc(fit FitFuncType) {
	pop.fit = fit
}

func (pop *DNAPopulation) SetConfig(size, varianceDim, objectiveDim int, lb, ub float64) {
	pop.size = size
	pop.varianceDim = varianceDim
	pop.objectiveDim = objectiveDim
	pop.lb = lb
	pop.ub = ub
}

func (pop *DNAPopulation) Size() int {
	return len(pop.individuals)
}

func (pop *DNAPopulation) Init() {
	agents := make([]*DNAAgent, pop.size)
	for i := range agents {
		agents[i] = CreateDNAAgent(pop.varianceDim, pop.lb, pop.ub)
	}
	pop.individuals = agents
}

func (pop *DNAPopulation) Append(invs []*DNAAgent) {
	if len(invs) > 0 {
		pop.individuals = append(pop.individuals, invs...)
		pop.size = len(pop.individuals)
	}
}

func (pop *DNAPopulation) Fit() [][]float64 {
	fits, ZMin := pop.fit(pop.individuals)
	pop.zmin = ZMin
	return fits
}

func (pop *DNAPopulation) ZMin() []float64 {
	return pop.zmin
}

func (pop *DNAPopulation) At(i int) *DNAAgent {
	return pop.individuals[i]
}

func (pop *DNAPopulation) VarianceDim() int {
	return pop.varianceDim
}

func (pop *DNAPopulation) ObjectiveDim() int {
	return pop.objectiveDim
}

func (pop *DNAPopulation) UpdatePosition(i int, position []float64) {
	for i, v := range position {
		position[i] = math.Round(v)
	}
	pop.individuals[i].UpdatePosition(position)
}

func (pop *DNAPopulation) LB() float64 {
	return pop.lb
}

func (pop *DNAPopulation) UB() float64 {
	return pop.ub
}

func (pop *DNAPopulation) PostWork() {
	for i := range pop.Size() {
		pop.individuals[i].PostWork()
	}
}
