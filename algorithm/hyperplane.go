package algorithm

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

const w float64 = 1e6

func asf(n, m int, localPoints [][]float64) []float64 {
	extrems := []float64{}
	for j := 0; j < m; j++ {
		minASFIndexDimM := 0
		minASF := -1.0
		for i := 0; i < n; i++ {
			maxASF := float64(0)
			for jj := 0; jj < m; jj++ {
				if j == jj {
					maxASF = max(maxASF, localPoints[i][jj])
				} else {
					maxASF = max(maxASF, localPoints[i][jj]*w)
				}

			}
			if minASFIndexDimM == 0 || maxASF < minASF {
				minASF = maxASF
				minASFIndexDimM = i
			}
		}
		extrems = append(extrems, localPoints[minASFIndexDimM]...)

	}
	//for i := 0; i < n; i++ {
	//	for j := 0; j < m; j++ {
	//		localPoints[i][j] -= ZMin[j]
	//	}
	//}

	one := make([]float64, m)
	for i := range one {
		one[i] = 1
	}
	ones := mat.NewDense(m, 1, one)

	extMtrix := mat.NewDense(m, m, extrems)
	var aInv mat.Dense
	aInv.Solve(extMtrix, ones)

	intercept := make([]float64, m)
	aInv.T()
	for j := 0; j < m; j++ {
		intercept[j] = aInv.At(j, 0)
	}
	for i := range intercept {
		if intercept[i] == 0 {
			for j := range localPoints {
				if intercept[i] < localPoints[j][i] {
					intercept[i] = localPoints[j][i]
				}
			}
		} else {
			intercept[i] = 1.0 / intercept[i]
		}
	}
	return intercept
}

// points: indicators matrix
// ZMin: global min value for each dimension
func bestIndex(points [][]float64, ZMin []float64, norm bool) (int, []float64, error) {
	n := len(points)
	if n == 0 {
		return 0, nil, errors.New("There is no points")
	}
	m := len(points[0])
	if m == 0 {
		return 0, nil, errors.New("Dimension is zero")
	}

	localPoints := make([][]float64, n)
	for i := range localPoints {
		copy(localPoints[i], points[i])
	}
	copy(localPoints, points)

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			localPoints[i][j] -= ZMin[j]
		}
	}

	// ASF
	intercept := asf(n, m, localPoints)

	// 点到平面距离，并不是真的距离，只是不影响比较大小
	distances := make([]float64, len(localPoints))
	minDistance := -1.1
	index := 0
	for i := range localPoints {
		distance := float64(0)
		for j := 0; j < m; j++ {
			if norm {
				if intercept[j] != 0 {
					distance += localPoints[i][j] / intercept[j]
				}
			} else {
				distance += intercept[j] * localPoints[i][j]
			}

		}
		distances[i] = distance
		if index == 0 || distance < minDistance {
			index = i
			minDistance = distance
		}
	}
	return index, distances, nil
}

func distanceToPlane(intercept []float64, point []float64) float64 {
	distance := float64(0)
	for j := range intercept {
		distance += intercept[j] * point[j]
	}
	return distance
}
