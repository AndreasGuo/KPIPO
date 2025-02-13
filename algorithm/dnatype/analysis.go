package DNAType

import (
	"fmt"
	"math"
	"slices"
	"strings"
)

func Continuity(seq Seq, t int) int {
	last := -1
	count := 1
	value := 0

	for _, base := range seq {
		if base == last {
			count += 1
		} else {
			if count >= t {
				value += count * count
			}
			count = 1
		}
		last = base
	}

	// continuous at tail
	if count >= t {
		value += count * count
	}

	return value
}

// define two bases are complementary or not
func cp(base1, base2 int) int {
	if base1+base2 == 3 {
		return 1
	}
	return 0
}

// threshold
func th(i, threshold int) int {
	if i > threshold {
		return i
	}
	return 0
}
func Hairpin(seq Seq, rMin int, pMin, pinLenT int) int {
	value := 0
	l := len(seq)
	// p is length of stem
	for p := pMin; p <= (l-rMin)/2; p++ {
		// r is length of loop
		for r := rMin; r <= l-2*p; r++ {
			// "i" is length of unused bases
			for i := 0; i <= l-2*p-r; i++ {
				sCp := 0
				// pin len defined shorter stem of a hairpin
				// suppose a real hair pin, it has a long side and a short side
				// which p+i is the long side
				pinLen := min(p+i, l-p-i-r)
				// check if complementary from both side of loop
				for j := 0; j < pinLen; j++ {
					sCp += cp(seq[p+i-j-1], seq[p+i+r+j])
				}
				value += th(sCp, pinLen/pinLenT)
			}
		}
	}
	return value
}

func shift(seq Seq, i int) Seq {
	// right shift i>0; left shift i<0
	l := len(seq)
	if i == 0 {
		return seq
	}
	absI := int(math.Abs(float64(i)))
	temp := make([]int, absI)
	for i := range temp {
		temp[i] = 5
	}
	switch {
	case 0 < i && i < l:
		return append(temp, seq[0:l-i]...)
	case i < 0 && i > -l:
		return append(seq[absI:l], temp...)
	default:
		temp = make([]int, l)
		for i := range temp {
			temp[i] = 5
		}
		return temp
	}
}

func lenNB(seq Seq) int {
	count := 0
	for _, base := range seq {
		if base <= G {
			count += 1
		}
	}
	return count
}

const HDIS = 0.17

func hDis(seq1, seq2 Seq) int {
	// this requires length of seq2 not less than seq1
	if len(seq2) < len(seq1) {
		panic("length of seq2 less than seq1 in hDis.")
	}
	cpS := 0
	l := len(seq1)
	for i := 0; i < l; i++ {
		cpS += cp(seq1[i], seq2[i])
	}
	threshold := HDIS * float64(lenNB(seq2)) / 2
	return th(cpS, int(threshold))
}

// continuous complementary start from i
func ccp(seq1, seq2 Seq, i int) int {
	c := 0
	for j := i; j < len(seq1); j++ {
		if seq1[j] == G-seq2[j] {
			c += 1
		} else {
			break
		}
	}
	return c
}

func hCon(seq1, seq2 Seq) int {
	HCON := 6
	value := 0
	for i := 0; i < len(seq1); i++ {
		value = value + th(ccp(seq1, seq2, i), HCON)
	}
	return value
}

func HMeasure(lSeq, rSeq Seq) int {
	l := len(lSeq)
	gap := int(math.Round(float64(l) / 4))
	HM := 0
	reversedLSeq := make(Seq, len(lSeq))
	copy(reversedLSeq, lSeq)
	slices.Reverse(reversedLSeq)
	for g := 0; g <= gap; g++ {
		temp := make([]int, g)
		for i := range temp {
			temp[i] = 5
		}
		tempSeq := append(rSeq, temp...)
		tempSeq = append(tempSeq, rSeq...)

		for i := -l + 1; i < l-1; i++ {
			shiftedSeq := shift(tempSeq, i)
			tempHM := hDis(reversedLSeq, shiftedSeq) + hCon(reversedLSeq, shiftedSeq)
			HM = max(HM, tempHM)
		}
	}
	return HM
}

const SDIS = 0.17

func sDis(seq1, seq2 Seq) int {
	sEq := 0
	for i := 0; i < len(seq1); i++ {
		if seq1[i] == seq2[i] {
			sEq += 1
		}
	}
	return th(sEq, int(SDIS*float64(len(seq2))/2))
}

func cEq(seq1, seq2 Seq, i int) int {
	c := 0
	for j := i; j < len(seq1); j++ {
		if seq1[j] == seq2[j] {
			c += 1
		} else {
			break
		}
	}
	return c
}

const SCON = 6

func sCon(seq1, seq2 Seq) int {
	value := 0
	for i := 0; i < len(seq1); i++ {
		value += th(cEq(seq1, seq2, i), SCON)
	}
	return value
}

func Similarity(lSeq, rSeq Seq) int {
	l := len(lSeq)
	gap := int(math.Round(float64(l) / 4))
	SM := 0
	//slices.Reverse(lSeq)
	for g := 0; g < gap; g++ {
		temp := make([]int, g)
		for i := range temp {
			temp[i] = 5
		}
		tempSeq := append(rSeq, temp...)
		tempSeq = append(tempSeq, rSeq...)

		for i := -l + 1; i < l; i++ {
			shiftedSeq := shift(tempSeq, i)
			tempSm := sDis(lSeq, shiftedSeq) + sCon(lSeq, shiftedSeq)
			SM = max(SM, tempSm)
		}
	}
	return SM
}

// 与matlab计算出的结果虽然不太一样但是可以接近
// 我觉得问题在于是此包使用的是
// Correction for deltaS: 0.368 x (N-1) x ln[Na+]
// (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)
// 问题不大，可以用，如果需要准确数值就再去matlab跑一下就好了
func MeltingTemperature(seq Seq) float64 {
	DNACode, _ := seq.ToStr()
	salt := 1.0
	primerConc := 1e-8
	return calcTm(DNACode, salt, primerConc)
}

var nearestNeighborParams = map[string]struct {
	deltaH float64 // 焓变（kcal/mol）
	deltaS float64 // 熵变（cal/mol/K）
}{
	"AA": {-7.9, -22.2},
	"AC": {-8.4, -22.4},
	"AG": {-7.8, -21},
	"AT": {-7.2, -20.4},
	"CA": {-8.5, -22.7},
	"CC": {-8.0, -19.9},
	"CG": {-10.6, -27.2},
	"CT": {-7.8, -21},
	"GA": {-8.2, -22.2},
	"GC": {-9.8, -24.4},
	"GG": {-8.0, -19.9},
	"GT": {-8.4, -22.4},
	"TA": {-7.2, -21.3},
	"TC": {-8.2, -22.2},
	"TG": {-8.5, -22.7},
	"TT": {-7.9, -22.2},
}

// 计算邻近配对法的熔解温度（Tm），并考虑引物浓度
func calcTm(sequence string, naConcentration, primerConcentration float64) float64 {
	var deltaH, deltaS float64
	seqLen := len(sequence)
	allSelfCompFlag := selfComp(sequence)
	// 处理每一对碱基
	for i := 0; i < seqLen-1; i++ {
		pair := strings.ToUpper(sequence[i : i+2]) // 取相邻的两个碱基对
		if params, exists := nearestNeighborParams[pair]; exists {
			deltaH += params.deltaH
			deltaS += params.deltaS
		} else {
			// 如果遇到不在参数表中的碱基对，跳过
			fmt.Printf("Warning: Invalid base pair %s\n", pair)
		}
	}

	if allSelfCompFlag {
		deltaS += -1.4
	}

	if sequence[0] == 'C' || sequence[0] == 'G' {
		deltaH += 0.1
		deltaS += -2.8
	} else if sequence[0] == 'A' || sequence[0] == 'T' {
		deltaH += 2.3
		deltaS += 4.1
	}

	if sequence[seqLen-1] == 'C' || sequence[seqLen-1] == 'G' {
		deltaH += 0.1
		deltaS += -2.8
	} else if sequence[seqLen-1] == 'A' || sequence[seqLen-1] == 'T' {
		deltaH += 2.3
		deltaS += 4.1
	}

	// 计算熔解温度 Tm
	// Tm = ΔH / (ΔS + R * ln[Na+])
	// R = 1.987 cal/mol/K
	R := 1.9872
	b := 4.0
	if allSelfCompFlag {
		b = 1.0
	}
	tm := deltaH * 1000 / (deltaS + R*math.Log(primerConcentration/b))
	tm += 16.0 * math.Log10(naConcentration)

	tm -= 273.15 // K to C
	return tm
}

func selfComp(sequence string) bool {
	seqLen := len(sequence)
	for i := range int(seqLen / 2) {
		if sequence[i] == 'A' && sequence[seqLen-1-i] != 'T' {
			return false
		}
		if sequence[i] == 'C' && sequence[seqLen-1-i] != 'G' {
			return false
		}
		if sequence[i] == 'G' && sequence[seqLen-1-i] != 'C' {
			return false
		}
		if sequence[i] == 'T' && sequence[seqLen-1-i] != 'A' {
			return false
		}
	}
	return true
}
