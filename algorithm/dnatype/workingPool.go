package DNAType

type SeqMapSingle struct {
	Index int
	Seq   Seq
}
type ResultMapSingle struct {
	Index int
	Value float64
}
type SeqMapPair struct {
	Index1, Index2 int
	Seq1, Seq2     Seq
}
type ResultMapPair struct {
	Index1, Index2 int
	Value          float64
}

func continuityWorker(in <-chan SeqMapSingle, out chan<- ResultMapSingle) {
	for job := range in {
		continuity := Continuity(job.Seq, 3)
		out <- ResultMapSingle{job.Index, float64(continuity)}
	}
}

func hairpinWorker(in <-chan SeqMapSingle, out chan<- ResultMapSingle) {
	for job := range in {
		hairpin := Hairpin(job.Seq, 6, 6, 3)
		out <- ResultMapSingle{job.Index, float64(hairpin)}
	}
}

func hmeasureWorker(in <-chan SeqMapPair, out chan<- ResultMapPair) {
	for job := range in {
		hm := HMeasure(job.Seq1, job.Seq2)
		out <- ResultMapPair{job.Index1, job.Index2, float64(hm)}
	}
}

func similarityWorker(in <-chan SeqMapPair, out chan<- ResultMapPair) {
	for job := range in {
		sm := Similarity(job.Seq1, job.Seq2)
		out <- ResultMapPair{job.Index1, job.Index2, float64(sm)}
	}
}

func meltingTemperatureWorker(in <-chan SeqMapSingle, out chan<- ResultMapSingle) {
	for job := range in {
		mt := MeltingTemperature(job.Seq)
		out <- ResultMapSingle{job.Index, float64(mt)}
	}
}

type FitChan struct {
	CtIn chan SeqMapSingle
	CtRe chan ResultMapSingle
	HpIn chan SeqMapSingle
	HpRe chan ResultMapSingle
	HmIn chan SeqMapPair
	HmRe chan ResultMapPair
	SmIn chan SeqMapPair
	SmRe chan ResultMapPair
	MtIn chan SeqMapSingle
	MtRe chan ResultMapSingle
}

func (fitChan *FitChan) Close() {
	close(fitChan.CtIn)
	close(fitChan.CtRe)
	close(fitChan.HpIn)
	close(fitChan.HpRe)
	close(fitChan.HmIn)
	close(fitChan.HmRe)
	close(fitChan.SmIn)
	close(fitChan.SmRe)
	close(fitChan.MtIn)
	close(fitChan.MtRe)
}

func CreateWorker(numOfSingle, numOfPair int, bufferSize int) *FitChan {
	continuityChan := make(chan SeqMapSingle, bufferSize)
	continuityResult := make(chan ResultMapSingle, bufferSize)
	hairpinChan := make(chan SeqMapSingle, bufferSize)
	hairpinResult := make(chan ResultMapSingle, bufferSize)
	hmChan := make(chan SeqMapPair, bufferSize)
	hmResult := make(chan ResultMapPair, bufferSize)
	smChan := make(chan SeqMapPair, bufferSize)
	smResult := make(chan ResultMapPair, bufferSize)
	mtChan := make(chan SeqMapSingle, bufferSize)
	mtResult := make(chan ResultMapSingle, bufferSize)
	for i := 0; i < numOfSingle; i++ {
		go continuityWorker(continuityChan, continuityResult)
		go hairpinWorker(hairpinChan, hairpinResult)
		go meltingTemperatureWorker(mtChan, mtResult)
	}
	for i := 0; i < numOfPair; i++ {
		go hmeasureWorker(hmChan, hmResult)
		go similarityWorker(smChan, smResult)
	}
	return &FitChan{continuityChan,
		continuityResult,
		hairpinChan,
		hairpinResult,
		hmChan,
		hmResult,
		smChan,
		smResult,
		mtChan,
		mtResult}
}
