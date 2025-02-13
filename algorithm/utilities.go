package algorithm

import (
	"math"
	"math/rand"
)

func sum[T int | float64](in []T) T {
	var s T = 0
	for i := range in {
		s += in[i]
	}
	return s
}

// 柯西分布
func C(gamma, x0 float64) float64 {
	r := rand.Float64()*0.9 + 0.05
	y := r - 0.5
	y *= math.Pi
	y = math.Tan(y)*gamma + x0
	return y
}

func C0(gamma float64) float64 {
	return C(gamma, 0)
}
