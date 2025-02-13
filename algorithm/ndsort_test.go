package algorithm

import "testing"

func TestNDSort(t *testing.T) {
	data := [][]int{
		{0, 1, 1, 1},
		{0, 0, 5, 3},
		{0, 0, 6, 3},
	}
	t.Log(NDSort(data))
}
