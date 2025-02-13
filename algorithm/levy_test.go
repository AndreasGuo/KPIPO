package algorithm

import (
	"fmt"
	"testing"
)

func TestLevy(t *testing.T) {
	l := levy(20)
	for i := range l {
		fmt.Println(l[i])
	}
}
