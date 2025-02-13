package DNAType

import (
	"slices"
	"sync"
)

type DNASet []DNAAgent
type FitFuncType func([]*DNAAgent) ([][]float64, []float64)

// 注意fit所给mt已经是偏离度了
func FitnessCall(dnaSet []*DNAAgent, index int, fitChan *FitChan, minVal float64, reversed bool) FitFuncType {
	var minCt, minHp, minHm, minSm, minMT = 400.0, 400.0, 400.0, 400.0, 400.0
	seqSet := make([]Seq, len(dnaSet))
	for i := range dnaSet {
		represent := dnaSet[i].Represent()
		seqSet[i] = make(Seq, len(represent))
		copy(seqSet[i], represent)
	}

	var mtValues = make([]float64, len(dnaSet))
	for i, seq := range seqSet {
		fitChan.MtIn <- SeqMapSingle{Index: i, Seq: seq}
		re := <-fitChan.MtRe
		mtValues[re.Index] = re.Value
	}

	return func(invs []*DNAAgent) ([][]float64, []float64) {
		go func() {
			for i := range invs {
				fitChan.CtIn <- SeqMapSingle{Index: i, Seq: invs[i].Represent()}
			}
		}()
		go func() {
			for i := range invs {
				fitChan.HpIn <- SeqMapSingle{Index: i, Seq: invs[i].Represent()}
			}
		}()
		go func() {
			for i := range invs {
				for j := range seqSet {
					if j == index {
						fitChan.HmIn <- SeqMapPair{Index1: i, Index2: j, Seq1: invs[i].Represent(), Seq2: invs[i].Represent()}
					} else {
						if !reversed {
							// 正常的算法
							fitChan.HmIn <- SeqMapPair{Index1: i, Index2: j, Seq1: invs[i].Represent(), Seq2: seqSet[j]}
						} else {
							// 交换前后
							fitChan.HmIn <- SeqMapPair{Index1: i, Index2: j, Seq1: seqSet[j], Seq2: invs[i].Represent()}
						}
					}
				}
			}
		}()
		go func() {
			for i := range invs {
				for j := range seqSet {
					if j == index {
						continue
					} else {
						if !reversed {
							// 正常的算法
							fitChan.SmIn <- SeqMapPair{Index1: i, Index2: j, Seq1: invs[i].Represent(), Seq2: seqSet[j]}
						} else {
							// 交换前后
							fitChan.SmIn <- SeqMapPair{Index1: i, Index2: j, Seq1: seqSet[j], Seq2: invs[i].Represent()}
						}
					}
				}
			}
		}()
		go func() {
			for i := range invs {
				fitChan.MtIn <- SeqMapSingle{Index: i, Seq: invs[i].Represent()}
			}
		}()

		continuityList := make([]float64, len(invs))
		hairpinList := make([]float64, len(invs))
		hmTable := make([][]float64, len(invs))
		for i := range hmTable {
			hmTable[i] = make([]float64, len(dnaSet))
		}
		smTable := make([][]float64, len(invs))
		for i := range smTable {
			smTable[i] = make([]float64, len(dnaSet))
		}
		mtList := make([]float64, len(invs))

		// 等鸡啄完了米山、狗舔完了面山，火烧断了金锁
		group := sync.WaitGroup{}
		group.Add(5)
		go func() {
			for i := 0; i < len(invs); i++ {
				ctRe := <-fitChan.CtRe
				continuityList[ctRe.Index] = ctRe.Value
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				hpRe := <-fitChan.HpRe
				hairpinList[hpRe.Index] = hpRe.Value
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				for j := 0; j < len(seqSet); j++ {
					hmRe := <-fitChan.HmRe
					hmTable[hmRe.Index1][hmRe.Index2] = hmRe.Value
				}
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				for j := 0; j < len(seqSet); j++ {
					if j == index {
						continue
					}
					smRe := <-fitChan.SmRe
					smTable[smRe.Index1][smRe.Index2] = smRe.Value
				}
			}
			group.Done()
		}()
		go func() {
			for i := 0; i < len(invs); i++ {
				mtRe := <-fitChan.MtRe
				mtList[mtRe.Index] = mtRe.Value
			}
			group.Done()
		}()
		group.Wait()

		hmList := sumTableRow(hmTable)
		smList := sumTableRow(smTable)
		// 连续与发卡使用实际值即可
		// 但是对于H-测度与相似性，一条DNA链的改变会影响其他链的值
		// 考虑使用
		// 1、七条的总值 || 平均值
		// 2、它对于其他DNA链在两个数值上的贡献度
		minCt = min(minCt, slices.Min(continuityList))
		minHp = min(minHp, slices.Min(hairpinList))
		minHm = min(minHm, slices.Min(hmList))
		minSm = min(minSm, slices.Min(smList))

		mtDiviant(mtList, mtValues, minVal)
		minMT = min(minMT, slices.Min(mtList))

		fits := [][]float64{}
		for i := 0; i < len(invs); i++ {
			agentFit := []float64{continuityList[i], hairpinList[i], hmList[i], smList[i], mtList[i]}
			invs[i].SetObjs(agentFit)
			fits = append(fits, agentFit)
		}
		return fits, []float64{minCt, minHp, minHm, minSm, minMT}
	}
}

func mtDiviant(values, compared []float64, minVal float64) {
	var sumCompared float64 = 0
	for _, value := range compared {
		sumCompared += value
	}
	avg := sumCompared / float64(len(compared))
	for i, value := range values {
		//avg := (sumCompared + value) / float64(len(compared)+1)
		//var div float64 = 0
		//minVal := min(value, slices.Min(compared))
		//maxVal := max(value, slices.Max(compared))
		//for _, c := range compared {
		//	div += math.Pow((c-avg)/(maxVal-minVal), 2)
		//}
		//div += math.Pow((value-avg)/(maxVal-minVal), 2)
		//div /= float64(config.DNASIZE)
		values[i] = max(value-avg, 0.5) // minVal)
	}
}

func sumTableRow[T int | float64](table [][]T) []float64 {
	sumList := make([]float64, len(table))
	for i, row := range table {
		var sum T = 0
		for _, element := range row {
			sum += element
		}
		sumList[i] = float64(sum)
	}
	return sumList
}

func avgTableRow[T int | float64](table [][]T) []float64 {
	avgs := make([]float64, len(table))
	for i, row := range table {
		var sum T = 0
		for _, element := range row {
			sum += element
		}
		avgs[i] = float64(sum) / float64(len(row))
	}
	return avgs
}

func norm(values []float64, minVal, maxVal float64, MINVALUE float64) {
	for i := range values {
		values[i] = (values[i] - minVal) / (maxVal - minVal)
		values[i] = max(values[i], MINVALUE)
	}
}

func sum(l []float64) float64 {
	var s float64
	for i := range l {
		s += l[i]
	}
	return s
}
