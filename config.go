package main

type Config struct {
	DIM             int
	DNASIZE         int
	POPSIZE         int
	MAXIT           int
	LB              int
	UB              int
	DNASETITERATION int
	MINVALUE        float64
	FITREVERSE      bool
	PLANENORM       bool
	CHOOSETOOPT     int
	// ORIGINPO        bool
	CD bool
}

// func DefaultConfig() *Config {
// 	config := Config{20, 7, 50, 200, 0, 3, 1000, 2e-2}
// 	return &config
// }
