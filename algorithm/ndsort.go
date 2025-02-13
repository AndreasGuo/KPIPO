package algorithm

import (
	"errors"
	"slices"

	"github.com/mkmik/argsort"
)

const inf float64 = 0xfff

// NDSort
// Efficient non-dominated sort with sequential search (ENS-SS)
// @Input objs matrix of object values, size*object
// @Output front no of each individual
func NDSort[T int | float64](objs [][]T) []int {
	var maxFrontNo = -1
	var frontNo = make([]int, len(objs))
	for i := range frontNo {
		frontNo[i] = -1
	}

	for notInitialValue(frontNo) < len(objs) {
		maxFrontNo += 1
		//marked := false
		for i := range objs {
			if frontNo[i] == -1 {
				dominated := false

				for j := 0; j < len(objs); j++ {
					if i == j {
						continue
					}
					if (frontNo[j] == -1 || frontNo[j] == maxFrontNo) && dominate(objs[j], objs[i]) {
						dominated = true
						break
					}
					//if !marked && frontNo[j] == -1 && dominate(objs[j], objs[i]) {
					//	dominated = true
					//	break
					//} else if (frontNo[j] == maxFrontNo) && dominate(objs[j], objs[i]) {
					//	dominated = true
					//	break
					//}
				}

				if !dominated {
					frontNo[i] = maxFrontNo
					//marked = true
				}
			}
		}
	}
	return frontNo
}

func notInitialValue(frontNo []int) int {
	count := 0
	for _, v := range frontNo {
		if v != -1 {
			count += 1
		}
	}
	return count
}

func dominate[T int | float64](lhs, rhs []T) bool {
	count := 0
	atLeastLess := false
	for i := range lhs {
		if lhs[i] <= rhs[i] {
			count += 1
		}
		if lhs[i] < rhs[i] {
			atLeastLess = true
		}
	}
	return (count == len(lhs)) && atLeastLess
}

func extractMatrix(fits [][]float64, positions []int) [][]float64 {
	out := make([][]float64, len(positions))
	for i := range positions {
		length := len(fits[positions[i]])
		out[i] = make([]float64, length)
		copy(out[i], fits[positions[i]])
	}
	return out
}

// returns the best and all selected (index)
func NDKPSort(fits [][]float64, ZMin []float64, size int, norm bool, cd bool) (int, []int) {
	NDFront := NDSort(fits)
	maxFront := slices.Max(NDFront) + 1
	var zeroPosition []int
	selectedPosition := make([]int, 0)
	maxSelectedFront := -1
	for f := range maxFront {
		front := findPosition(NDFront, f)
		if f == 0 {
			zeroPosition = make([]int, len(front))
			copy(zeroPosition, front)
		}
		if len(selectedPosition)+len(front) < size {
			selectedPosition = append(selectedPosition, front...)
			maxSelectedFront = f
		} else if len(selectedPosition)+len(front) == size {
			selectedPosition = append(selectedPosition, front...)
			maxSelectedFront = f
			break
		} else {
			lastFrontSelected := selectFromLastFront(fits, front, ZMin, size-len(selectedPosition), norm)
			selectedPosition = append(selectedPosition, lastFrontSelected...)
			maxSelectedFront = f
		}
	}

	if maxSelectedFront == -1 {
		panic("??")
	}

	if maxSelectedFront == 0 {
		return selectedPosition[0], selectedPosition
	}

	zeroFront := extractMatrix(fits, zeroPosition)
	var bI int
	if !cd {
		bI, _, _ = bestIndex(zeroFront, ZMin, norm)
	} else {
		normedObjs := normObjs(zeroFront)
		cds, _ := crowdingDistance(normedObjs)
		bI = maxCrowdingDistanceIndex(cds)
	}

	return zeroPosition[bI], selectedPosition
}

func selectFromLastFront(fits [][]float64, front []int, ZMin []float64, remainSize int, norm bool) []int {
	frontFits := make([][]float64, len(front))
	for i := range front {
		frontFits[i] = make([]float64, len(fits[0]))
		copy(frontFits[i], fits[front[i]])
	}
	_, distances, _ := bestIndex(frontFits, ZMin, norm)
	indicies := argsort.SortSlice(distances, func(i, j int) bool {
		return distances[i] < distances[j]
	})
	selected := make([]int, 0)
	for i := range remainSize {
		selected = append(selected, front[indicies[i]])
	}
	return selected
}

func findPosition(lhs []int, rhs int) []int {
	positions := []int{}
	for i := range lhs {
		if lhs[i] == rhs {
			positions = append(positions, i)
		}
	}
	return positions
}

func normObjs(fits [][]float64) [][]float64 {
	n := len(fits)
	normedFits := make([][]float64, n)
	if n == 0 {
		return normedFits
	}
	m := len(fits[0])
	for i := range n {
		normedFits[i] = make([]float64, m)
		interval := max(slices.Max(fits[i])-slices.Min(fits[i]), 1e-2)
		for j := range m {
			normedFits[i][j] = fits[i][j] / interval
		}
	}

	return normedFits
}

func crowdingDistance(normedFits [][]float64) ([]float64, error) {
	// each column
	n := len(normedFits)
	if n == 0 {
		return nil, errors.New("no fits values")
	}
	m := len(normedFits[0])
	cd := make([]float64, n)
	for j := range m {
		obj_m := make([]float64, 0, n)
		for i := range n {
			obj_m = append(obj_m, normedFits[i][j])
		}
		indicies := argsort.SortSlice(obj_m, func(i, j int) bool {
			return obj_m[i] < obj_m[j]
		})
		interval := max(slices.Max(obj_m)-slices.Min(obj_m), 0.5)
		for i := range n {
			if indicies[i] == 0 || indicies[i] == n-1 {
				cd[i] += inf
			} else {
				cd[i] += (obj_m[indicies[i]-1] - obj_m[indicies[i]+1]) / interval
			}
		}
	}
	return cd, nil
}

func maxCrowdingDistanceIndex(crowdingDistance []float64) (index int) {
	index = 0
	for i := range len(crowdingDistance) {
		if crowdingDistance[index] > inf {
			break
		}
		if crowdingDistance[i] > crowdingDistance[index] {
			index = i
		}
	}
	return
}
