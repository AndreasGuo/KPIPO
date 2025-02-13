package DNAType

import (
	"errors"
	"strings"
)

const C = 0
const T = 1
const A = 2
const G = 3

type DNASeq interface {
	ToStr() (string, error)
}

type Seq []int

func (seq Seq) ToStr() (string, error) {
	builder := strings.Builder{}
	builder.Grow(len(seq))
	for _, base := range seq {
		switch base {
		case C:
			builder.WriteByte('C')
		case A:
			builder.WriteByte('A')
		case T:
			builder.WriteByte('T')
		case G:
			builder.WriteByte('G')
		default:
			return "", errors.New("undefined base token")
		}
	}
	return builder.String(), nil
}

func ToSeq(str string) (Seq, error) {
	seq := make(Seq, len(str))
	for i := range str {
		switch str[i] {
		case 'C':
			seq[i] = C
		case 'A':
			seq[i] = A
		case 'T':
			seq[i] = T
		case 'G':
			seq[i] = G
		default:
			return seq, errors.New("undefined base name")
		}
	}
	return seq, nil
}
