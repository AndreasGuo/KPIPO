package DNAType

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

type DNAAgent struct {
	DNARowVariance []float64
	Seq            Seq
	Ct             float64
	Hp             float64
	Hm             float64
	Sm             float64
	Mt             float64
}

func (dnaAgent *DNAAgent) SetObjs(objs []float64) {
	dnaAgent.Ct = objs[0]
	dnaAgent.Hp = objs[1]
	dnaAgent.Hm = objs[2]
	dnaAgent.Sm = objs[3]
	dnaAgent.Mt = objs[4]
}

func (dnaAgent *DNAAgent) Clone() *DNAAgent {
	variance := dnaAgent.DNARowVariance
	seq := dnaAgent.Seq

	newVariance := make([]float64, len(variance))
	newSeq := make(Seq, len(seq))

	copy(newVariance, variance)
	copy(newSeq, seq)
	agent := DNAAgent{newVariance,
		newSeq,
		dnaAgent.Ct,
		dnaAgent.Hp,
		dnaAgent.Hm,
		dnaAgent.Sm,
		dnaAgent.Mt,
	}
	return &agent
}

func (dnaAgent *DNAAgent) String() (string, error) {
	return dnaAgent.Seq.ToStr()
}

func (dnaAgent *DNAAgent) Represent() []int {
	return dnaAgent.Seq
}

func (dnaAgent *DNAAgent) UpdatePosition(position []float64) {
	dnaAgent.DNARowVariance = position
}

func (dnaAgent *DNAAgent) Objs() []float64 {
	objs := make([]float64, 5)
	objs[0] = float64(dnaAgent.Ct)
	objs[1] = float64(dnaAgent.Hp)
	objs[2] = float64(dnaAgent.Hm)
	objs[3] = float64(dnaAgent.Sm)
	objs[4] = float64(dnaAgent.Mt)
	return objs
}

func (dnaAgent *DNAAgent) Variance() []float64 {
	return dnaAgent.DNARowVariance
}

func CreateDNAAgent(dim int, lb, ub float64) *DNAAgent {
	agent := DNAAgent{}
	// init raw variance
	agent.DNARowVariance = make([]float64, dim)
	for i := 0; i < dim; i++ {
		agent.DNARowVariance[i] = math.Round(rand.Float64()*(ub-lb) + lb)
	}
	agent.RepairAndToSeq()
	return &agent
}

func (dnaAgent *DNAAgent) PostWork() {
	dnaAgent.RepairAndToSeq()
}

func (dnaAgent *DNAAgent) RepairAndToSeq() {
	dnaAgent.fixGCContent()
	//dnaAgent.NoRunLength()
	// generate DNA seq.
	dnaAgent.Seq = convSeq(dnaAgent.DNARowVariance)
}

func (dnaAgent *DNAAgent) fixGCContent() {
	// First control GC at 50%
	GCPosition := []int{} //make([]int, len(agent.variance))
	ATPosition := []int{} //make([]int, len(agent.variance))
	for i, value := range dnaAgent.DNARowVariance {
		if math.IsNaN(value) {
			dnaAgent.DNARowVariance[i] = 3
			value = 3
		}
		if value == C || value == G {
			GCPosition = append(GCPosition, i)
		} else {
			ATPosition = append(ATPosition, i)
		}
	}
	if len(GCPosition) > len(dnaAgent.DNARowVariance)/2 {
		toRepairLen := len(GCPosition) - len(dnaAgent.DNARowVariance)/2
		shuffle(GCPosition)
		toRepair := GCPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.DNARowVariance[toRepair[i]] = T
			} else {
				dnaAgent.DNARowVariance[toRepair[i]] = A
			}
		}
	}

	if len(ATPosition) > len(dnaAgent.DNARowVariance)/2 {
		toRepairLen := len(ATPosition) - len(dnaAgent.DNARowVariance)/2
		shuffle(ATPosition)
		toRepair := ATPosition[0:toRepairLen]
		for i := 0; i < toRepairLen; i++ {
			if rand.Float64() < 0.5 {
				dnaAgent.DNARowVariance[toRepair[i]] = C
			} else {
				dnaAgent.DNARowVariance[toRepair[i]] = G
			}
		}
	}
}

func (agent *DNAAgent) NoRunLength() {
	// Assure no continuous identical bases
	for i := 0; i < len(agent.DNARowVariance)-1; i++ {
		if agent.DNARowVariance[i] == agent.DNARowVariance[i+1] {
			// i+1位会被换掉
			for j := 0; j < len(agent.DNARowVariance); j++ {
				if i == j || i+1 == j {
					continue
				}
				// 要满足的条件：
				// 1 被换过来的碱基(j)与i不同
				// 2 j<i时要左右都不同，j>i时只保证与左不同
				if agent.DNARowVariance[j] != agent.DNARowVariance[i] {
					if j < i {
						if j == 0 && agent.DNARowVariance[i+1] != agent.DNARowVariance[j+1] {
							agent.DNARowVariance[j], agent.DNARowVariance[i+1] = agent.DNARowVariance[i+1], agent.DNARowVariance[j]
							break
						} else if agent.DNARowVariance[i+1] != agent.DNARowVariance[j+1] && agent.DNARowVariance[i+1] != agent.DNARowVariance[j-1] {
							agent.DNARowVariance[j], agent.DNARowVariance[i+1] = agent.DNARowVariance[i+1], agent.DNARowVariance[j]
							break
						}
					} else if j == i+2 {
						if agent.DNARowVariance[j] != agent.DNARowVariance[i+1] {
							agent.DNARowVariance[j], agent.DNARowVariance[i+1] = agent.DNARowVariance[i+1], agent.DNARowVariance[j]
							break
						}
					} else {
						if agent.DNARowVariance[i+1] != agent.DNARowVariance[j-1] {
							agent.DNARowVariance[j], agent.DNARowVariance[i+1] = agent.DNARowVariance[i+1], agent.DNARowVariance[j]
							break
						}
					}
				}
			}
		}
	}
}

func convSeq(variance []float64) Seq {
	// convert float variance array to DNA Seq
	seq := make(Seq, len(variance))
	for i, value := range variance {
		switch int(value) {
		case 0:
			seq[i] = C
		case 1:
			seq[i] = T
		case 2:
			seq[i] = A
		case 3:
			seq[i] = G
		default:
			panic(fmt.Sprintf("error while repair DNA seq, value %f do not defined!", value))
		}
	}
	return seq
}

func shuffle(slice []int) {
	r := rand.New(rand.NewSource(time.Now().Unix()))
	for len(slice) > 0 {
		n := len(slice)
		randIndex := r.Intn(n)
		slice[n-1], slice[randIndex] = slice[randIndex], slice[n-1]
		slice = slice[:n-1]
	}
}
