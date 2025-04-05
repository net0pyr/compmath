package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return 3*math.Pow(x, 3) + 1.7*math.Pow(x, 2) - 15.42*x + 6.89
}

func df(x float64) float64 {
	return 9*math.Pow(x, 2) + 3.4*x - 15.42
}

func bisectionMethod(a, b, epsilon float64) (float64, int) {
	var c float64
	iterations := 0
	for (b - a) >= epsilon {
		c = (a + b) / 2
		if f(c) == 0.0 {
			break
		} else if f(c)*f(a) < 0 {
			b = c
		} else {
			a = c
		}
		iterations++
	}
	return c, iterations
}

func chordMethod(a, b, epsilon float64) (float64, [][]float64) {
	var x float64
	table := [][]float64{}
	for {
		x = a - (f(a)*(b-a))/(f(b)-f(a))
		table = append(table, []float64{a, b, x, f(a), f(b), f(x), math.Abs(x - a)})
		if math.Abs(f(x)) < epsilon {
			break
		}
		if f(a)*f(x) < 0 {
			b = x
		} else {
			a = x
		}
	}
	return x, table
}

func simpleIteration(x0, epsilon float64) (float64, [][]float64) {
	phi := func(x float64) float64 {
		return -0.0396*math.Pow(x, 3) - 0.02244*math.Pow(x, 2) + 1.203544*x - 0.090948
	}

	var x1 float64
	table := [][]float64{}
	for {
		x1 = phi(x0)
		table = append(table, []float64{x0, x1, f(x1), math.Abs(x1 - x0)})
		if math.Abs(x1-x0) < epsilon {
			break
		}
		x0 = x1
	}
	return x1, table
}

func secantMethod(x0, x1, epsilon float64) (float64, [][]float64) {
	var x2 float64
	table := [][]float64{}
	for {
		x2 = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0))
		table = append(table, []float64{x0, x1, x2, f(x2), math.Abs(x2 - x1)})
		if math.Abs(x2-x1) < epsilon {
			break
		}
		x0 = x1
		x1 = x2
	}
	return x2, table
}

func newtonMethod(x0, epsilon float64) float64 {
	var x1 float64
	for {
		fx := f(x0)
		dfx := df(x0)
		x1 = x0 - fx/dfx
		if math.Abs(x1-x0) < epsilon {
			break
		}
		x0 = x1
	}
	return x1
}

func generateLatexChordTable(table [][]float64) string {
	latex := "\\begin{longtable}{|c|c|c|c|c|c|c|c|}\n"
	latex += "\\caption{Метод хорд}\n"
	latex += "\\hline\n"
	latex += "№ шага & $a$ & $b$ & $x$ & $f(a)$ & $f(b)$ & $f(x)$ & $x_{k+1} - x_k$ \\\\\n"
	latex += "\\hline\n"
	for i, row := range table {
		latex += fmt.Sprintf("%d & %.3f & %.3f & %.3f & %.6f & %.6f & %.6f & %.6f \\\\\n", i+1, row[0], row[1], row[2], row[3], row[4], row[5], row[6])
	}
	latex += "\\hline\n"
	latex += "\\end{longtable}"
	return latex
}

func generateLatexSimpleIterationTable(table [][]float64) string {
	latex := "\\begin{longtable}{|c|c|c|c|c|}\n"
	latex += "\\caption{Метод простой итерации}\n"
	latex += "\\hline\n"
	latex += "№ итерации & $x_k$ & $x_{k+1}$ & $f(x_{k+1})$ & $x_{k+1} - x_k$ \\\\\n"
	latex += "\\hline\n"
	for i, row := range table {
		latex += fmt.Sprintf("%d & %.6f & %.6f & %.6f & %.6f \\\\\n", i+1, row[0], row[1], row[2], row[3])
	}
	latex += "\\hline\n"
	latex += "\\end{longtable}"
	return latex
}

func generateLatexSecantTable(table [][]float64) string {
	latex := "\\begin{longtable}{|c|c|c|c|c|c|}\n"
	latex += "\\caption{Метод секущих}\n"
	latex += "\\hline\n"
	latex += "№ итерации & $x_{k-1}$ & $x_k$ & $x_{k+1}$ & $f(x_{k+1})$ & $x_{k+1} - x_k$ \\\\\n"
	latex += "\\hline\n"
	for i, row := range table {
		latex += fmt.Sprintf("%d & %.6f & %.6f & %.6f & %.6f & %.6f \\\\\n", i+1, row[0], row[1], row[2], row[3], row[4])
	}
	latex += "\\hline\n"
	latex += "\\end{longtable}"
	return latex
}

func main() {
	epsilon := 0.01

	aBisection := 0.0
	bBisection := 1.0
	root1, iterationsBisection := bisectionMethod(aBisection, bBisection, epsilon)
	fmt.Printf("Метод половинного деления: корень = %.3f, итераций = %d\n", root1, iterationsBisection)

	x0Newton := 1.0
	root2 := newtonMethod(x0Newton, epsilon)
	fmt.Printf("Метод Ньютона: %.3f\n", root2)

	x0 := 2.0
	root3, tableSimpleIter := simpleIteration(x0, epsilon)
	fmt.Printf("Метод простой итерации: %.3f\n", root3)

	x0 = -3.0
	_, tableSimpleIter = simpleIteration(x0, epsilon)
	a := 0.0
	b := 1.0
	_, tableChord := chordMethod(a, b, epsilon)
	x0Secant := 1.0
	x1Secant := 2.0
	_, tableSecant := secantMethod(x0Secant, x1Secant, epsilon)
	fmt.Printf("Метод секущих: %.3f\n", root3)
	fmt.Println(generateLatexSimpleIterationTable(tableSimpleIter))
	fmt.Println(generateLatexChordTable(tableChord))
	fmt.Println(generateLatexSecantTable(tableSecant))
}
