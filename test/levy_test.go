package test

import (
	"fmt"
	"testing"
)

func TestLevy(t *testing.T) {
	l := Levy(20)
	for range 999 {
		ll := Levy(20)
		for i := range 20 {
			l[i] += ll[i]
		}
	}
	for i := range 20 {
		l[i] /= 1000
	}
	fmt.Println(l)
}
