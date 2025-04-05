package main

import (
	"fmt"
	"image/color"
	"math"
	"os"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func f1(x float64) float64 {
	return math.Pow(x, 3) - 2*x - 5
}

func f2(x float64) float64 {
	return math.Sin(x) - x/2
}

func f3(x float64) float64 {
	return math.Exp(x) - 3*x
}

func f4(x float64) float64 {
	return math.Log(x) - 1
}

func f5(x float64) float64 {
	return math.Pow(x, 2) - 2
}

func system1(x, y float64) (float64, float64) {
	return math.Sin(x) + y - 1.2, 2*x + math.Cos(y) - 2
}

func system2(x, y float64) (float64, float64) {
	return x*x + y*y - 4, math.Exp(x) + y - 1
}

func system3(x, y float64) (float64, float64) {
	return x + y*y - 4, x*x + y - 3
}

func system4(x, y float64) (float64, float64) {
	return math.Sin(x+y) - 1.2*x, x*x + y*y - 1
}

func system5(x, y float64) (float64, float64) {
	return x*x - y - 1, x + y*y - 1
}

func bisectionMethod(f func(float64) float64, a, b, epsilon float64) (float64, int) {
	var c float64
	iterations := 0
	for math.Abs(b-a) > epsilon {
		c = (a + b) / 2
		if math.Abs(f(c)) < epsilon {
			break
		} else if f(c)*f(a) < 0 {
			b = c
		} else {
			a = c
		}
		iterations++
	}
	return (a + b) / 2, iterations
}

func newtonMethod(f func(float64) float64, df func(float64) float64, x0, epsilon float64) (float64, int) {
	x := x0
	iterations := 0
	for math.Abs(f(x)) > epsilon {
		x = x - f(x)/df(x)
		iterations++
	}
	return x, iterations
}

func simpleIterationMethod(f func(float64) float64, df func(float64) float64, a, b, x0, epsilon float64) (float64, int) {
	maxDeriv := findMaxDerivative(df, a, b)
	var lambda float64

	if maxDeriv > 0 {
		lambda = -1 / maxDeriv
	} else {
		lambda = 1 / math.Abs(maxDeriv)
	}

	phi := func(x float64) float64 {
		return x + lambda*f(x)
	}

	const h = 1e-8

	dphidx := (phi(x0+h) - phi(x0)) / h

	if math.Abs(dphidx) >= 1 {
		fmt.Printf("Предупреждение: условие сходимости не выполняется (q = %.4f)\n", math.Abs(dphidx))
		fmt.Printf("Частные производные: dphi/dx=%.4f\n", dphidx)
	} else {
		fmt.Printf("Условие сходимости выполняется (q = %.4f < 1)\n", math.Abs(dphidx))
	}

	x := x0
	iterations := 0
	for math.Abs(f(x)) > epsilon {
		x = phi(x)
		iterations++
		if iterations > 1000 {
			panic("Превышено максимальное число итераций (1000)")
		}

		if x < a || x > b {
			panic("Итерационный процесс вышел за границы отрезка")
		}
	}
	return x, iterations
}

func findMaxDerivative(df func(float64) float64, a, b float64) float64 {
	const steps = 100
	step := (b - a) / steps
	max := math.Abs(df(a))

	for x := a; x <= b; x += step {
		current := math.Abs(df(x))
		if current > max {
			max = current
		}
	}
	return max
}

func simpleIterationSystem(phi1, phi2 func(float64, float64) float64,
	x0, y0, epsilon float64,
	maxIterations int) (float64, float64, int, []float64, bool) {

	const h = 1e-8

	dphi1dx := (phi1(x0+h, y0) - phi1(x0, y0)) / h
	dphi1dy := (phi1(x0, y0+h) - phi1(x0, y0)) / h
	dphi2dx := (phi2(x0+h, y0) - phi2(x0, y0)) / h
	dphi2dy := (phi2(x0, y0+h) - phi2(x0, y0)) / h

	q1 := math.Abs(dphi1dx) + math.Abs(dphi1dy)
	q2 := math.Abs(dphi2dx) + math.Abs(dphi2dy)
	maxQ := math.Max(q1, q2)

	if maxQ >= 1 {
		fmt.Printf("Предупреждение: условие сходимости не выполняется (q = %.4f)\n", maxQ)
		fmt.Printf("Частные производные: dphi1/dx=%.4f, dphi1/dy=%.4f, dphi2/dx=%.4f, dphi2/dy=%.4f\n",
			dphi1dx, dphi1dy, dphi2dx, dphi2dy)
	} else {
		fmt.Printf("Условие сходимости выполняется (q = %.4f < 1)\n", maxQ)
	}

	x, y := x0, y0
	iterations := 0
	errors := []float64{}
	converged := false

	for iterations < maxIterations {
		iterations++

		xNew := phi1(x, y)
		yNew := phi2(x, y)

		errorX := math.Abs(xNew - x)
		errorY := math.Abs(yNew - y)
		errors = []float64{errorX, errorY}

		if errorX < epsilon && errorY < epsilon {
			converged = true
			break
		}

		x, y = xNew, yNew
	}

	if converged {
		dphi1dx_final := (phi1(x+h, y) - phi1(x, y)) / h
		dphi1dy_final := (phi1(x, y+h) - phi1(x, y)) / h
		dphi2dx_final := (phi2(x+h, y) - phi2(x, y)) / h
		dphi2dy_final := (phi2(x, y+h) - phi2(x, y)) / h

		q1_final := math.Abs(dphi1dx_final) + math.Abs(dphi1dy_final)
		q2_final := math.Abs(dphi2dx_final) + math.Abs(dphi2dy_final)
		maxQ_final := math.Max(q1_final, q2_final)

		if maxQ_final >= 1 {
			fmt.Printf("Предупреждение: в решении условие сходимости не выполняется (q = %.4f)\n", maxQ_final)
		}
	}

	return x, y, iterations, errors, converged
}

func hasRoot(f func(float64) float64, a, b float64) bool {
	return f(a)*f(b) <= 0
}

func plotFunction(f func(float64) float64, a, b float64, filename string) error {
	p := plot.New()

	p.Title.Text = "График функции"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"

	p.Add(plotter.NewGrid())

	points := make(plotter.XYs, 0)
	for x := a; x <= b; x += (b - a) / 1000 {
		y := f(x)
		if !math.IsNaN(y) && !math.IsInf(y, 0) {
			points = append(points, plotter.XY{X: x, Y: y})
		}
	}

	line, err := plotter.NewLine(points)
	if err != nil {
		return err
	}
	line.Color = plotutil.Color(0)
	p.Add(line)

	roots := findRoots(f, a, b)
	rootPoints := make(plotter.XYs, len(roots))
	for i, x := range roots {
		rootPoints[i] = plotter.XY{X: x, Y: 0}
	}

	if len(rootPoints) > 0 {
		rootsScatter, err := plotter.NewScatter(rootPoints)
		if err != nil {
			return err
		}
		rootsScatter.GlyphStyle.Color = plotutil.Color(1)
		rootsScatter.GlyphStyle.Radius = vg.Points(4)
		p.Add(rootsScatter)
		p.Legend.Add("Корни", rootsScatter)
	}

	if a <= 0 && 0 <= b {
		y := f(0)
		if !math.IsNaN(y) && !math.IsInf(y, 0) {
			oyPoint, err := plotter.NewScatter(plotter.XYs{{X: 0, Y: y}})
			if err != nil {
				return err
			}
			oyPoint.GlyphStyle.Color = plotutil.Color(2)
			oyPoint.GlyphStyle.Radius = vg.Points(4)
			p.Add(oyPoint)
			p.Legend.Add("f(0)", oyPoint)
		}
	}

	p.Legend.Top = true
	p.Legend.Left = true
	p.Legend.Add("Функция", line)

	if err := p.Save(8*vg.Inch, 6*vg.Inch, filename); err != nil {
		return err
	}
	return nil
}

func findRoots(f func(float64) float64, a, b float64) []float64 {
	const steps = 1000
	var roots []float64
	step := (b - a) / steps

	for x := a; x <= b; x += step {
		if f(x)*f(x+step) <= 0 {
			root, _ := bisectionMethod(f, x, x+step, 1e-6)
			roots = append(roots, root)
		}
	}

	return roots
}

func plotSystem(f1, f2 func(float64, float64) float64, xmin, xmax, ymin, ymax, solutionX, solutionY float64, filename string) error {
	p := plot.New()

	p.Title.Text = "График системы уравнений"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	p.Add(plotter.NewGrid())

	grid1 := unitGrid{
		f:    f1,
		xmin: xmin,
		xmax: xmax,
		ymin: ymin,
		ymax: ymax,
		rows: 100,
		cols: 100,
	}

	contours1 := plotter.NewContour(grid1, []float64{0}, palette.Rainbow(10, 1, 1, 1, 1, 1))
	redLine, _ := plotter.NewLine(plotter.XYs{})
	redLine.Color = color.RGBA{R: 255, A: 255}
	contours1.LineStyles[0].Color = redLine.Color
	p.Add(contours1)

	grid2 := unitGrid{
		f:    f2,
		xmin: xmin,
		xmax: xmax,
		ymin: ymin,
		ymax: ymax,
		rows: 100,
		cols: 100,
	}

	contours2 := plotter.NewContour(grid2, []float64{0}, palette.Rainbow(10, 1, 1, 1, 1, 1))
	blueLine, _ := plotter.NewLine(plotter.XYs{})
	blueLine.Color = color.RGBA{B: 255, A: 255}
	contours2.LineStyles[0].Color = blueLine.Color
	p.Add(contours2)

	if !math.IsNaN(solutionX) && !math.IsNaN(solutionY) {
		solutionPoints := make(plotter.XYs, 1)
		solutionPoints[0].X = solutionX
		solutionPoints[0].Y = solutionY
		solPoint, err := plotter.NewScatter(solutionPoints)
		if err != nil {
			return fmt.Errorf("ошибка создания точки решения: %v", err)
		}
		solPoint.GlyphStyle.Color = color.RGBA{G: 255, A: 255}
		solPoint.GlyphStyle.Radius = vg.Points(6)
		p.Add(solPoint)
		p.Legend.Add("Решение", solPoint)
	}

	p.Legend.Add("Уравнение 1", redLine)
	p.Legend.Add("Уравнение 2", blueLine)
	p.Legend.Top = true

	if err := p.Save(8*vg.Inch, 6*vg.Inch, filename); err != nil {
		return fmt.Errorf("ошибка сохранения графика: %v", err)
	}

	return nil
}

type unitGrid struct {
	f                      func(x, y float64) float64
	xmin, xmax, ymin, ymax float64
	rows, cols             int
}

func (g unitGrid) Dims() (c, r int) { return g.cols, g.rows }
func (g unitGrid) X(c int) float64 {
	return g.xmin + (g.xmax-g.xmin)*float64(c)/float64(g.cols-1)
}
func (g unitGrid) Y(r int) float64 {
	return g.ymin + (g.ymax-g.ymin)*float64(r)/float64(g.rows-1)
}
func (g unitGrid) Z(c, r int) float64 {
	x := g.X(c)
	y := g.Y(r)
	return g.f(x, y)
}

func main() {
	var choice int
	fmt.Println("Выберите задачу:")
	fmt.Println("1. Решить одно уравнение")
	fmt.Println("2. Решить систему нелинейных уравнений")
	fmt.Scan(&choice)

	if choice == 1 {
		var eqChoice int
		fmt.Println("Выберите уравнение:")
		fmt.Println("1. x^3 - 2x - 5")
		fmt.Println("2. sin(x) - x/2")
		fmt.Println("3. exp(x) - 3x")
		fmt.Println("4. ln(x) - 1")
		fmt.Println("5. x^2 - 2")
		fmt.Scan(&eqChoice)

		var f func(float64) float64
		var df func(float64) float64
		var d2f func(float64) float64
		var defaultA, defaultB float64

		switch eqChoice {
		case 1:
			f = f1
			df = func(x float64) float64 { return 3*math.Pow(x, 2) - 2 }
			d2f = func(x float64) float64 { return 6 * x }
			defaultA, defaultB = 1.0, 3.0
		case 2:
			f = f2
			df = func(x float64) float64 { return math.Cos(x) - 0.5 }
			d2f = func(x float64) float64 { return -math.Sin(x) }
			defaultA, defaultB = 1.0, 3.0
		case 3:
			f = f3
			df = func(x float64) float64 { return math.Exp(x) - 3 }
			d2f = func(x float64) float64 { return math.Exp(x) }
			defaultA, defaultB = 0.0, 2.0
		case 4:
			f = f4
			df = func(x float64) float64 { return 1 / x }
			d2f = func(x float64) float64 { return -1 / (x * x) }
			defaultA, defaultB = 1.0, 3.0
		case 5:
			f = f5
			df = func(x float64) float64 { return 2 * x }
			d2f = func(x float64) float64 { return 2.0 }
			defaultA, defaultB = 1.0, 2.0
		default:
			fmt.Println("Неверный выбор уравнения.")
			return
		}

		err := plotFunction(f, defaultA-2, defaultB+2, "function_preview.png")
		if err != nil {
			fmt.Println("Ошибка при построении графика:", err)
		} else {
			fmt.Println("График функции сохранен в файл function_preview.png")
			fmt.Println("Пожалуйста, посмотрите график для выбора подходящего отрезка [a,b]")
		}

		var a, b, epsilon float64
		fmt.Printf("Введите границы интервала [a, b] (рекомендуемый диапазон около [%v, %v]):\n", defaultA, defaultB)
		fmt.Scan(&a, &b)

		fmt.Println("Введите погрешность вычисления (epsilon):")
		fmt.Scan(&epsilon)

		if !hasRoot(f, a, b) {
			fmt.Println("На интервале нет корня или их несколько.")
			return
		}

		file, err := os.Create("results.txt")
		if err != nil {
			fmt.Println("Ошибка при создании файла:", err)
			return
		}
		defer file.Close()

		rootBisection, iterBisection := bisectionMethod(f, a, b, epsilon)
		fmt.Printf("\nМетод половинного деления: корень(x0) = %v, итераций = %d\n", rootBisection, iterBisection)
		fmt.Fprintf(file, "Метод половинного деления: корень = %v, итераций = %d\n", rootBisection, iterBisection)
		fmt.Printf("Проверка метода половинного деления: f(x0)=%v\n\n", f(rootBisection))

		var x0 float64
		if f(a)*d2f(a) > 0 {
			x0 = a
		} else {
			x0 = b
		}
		fmt.Printf("Начальное приближение: %v\n\n", x0)

		rootNewton, iterNewton := newtonMethod(f, df, x0, epsilon)
		fmt.Printf("Метод Ньютона: корень(x0) = %v, итераций = %d\n", rootNewton, iterNewton)
		fmt.Fprintf(file, "Метод Ньютона: корень = %v, итераций = %d\n", rootNewton, iterNewton)
		fmt.Printf("Проверка метода Ньютона: f(x0)=%v\n\n", f(rootNewton))

		rootSimpleIter, iterSimpleIter := simpleIterationMethod(f, df, a, b, x0, epsilon)
		fmt.Printf("Метод простой итерации: корень = %v, итераций = %d\n", rootSimpleIter, iterSimpleIter)
		fmt.Fprintf(file, "Метод простой итерации: корень(x0) = %v, итераций = %d\n", rootSimpleIter, iterSimpleIter)
		fmt.Printf("Проверка простой итерации: f(x0)=%v\n\n", f(rootSimpleIter))

	} else if choice == 2 {
		var sysChoice int
		fmt.Println("Выберите систему уравнений:")
		fmt.Println("1. sin(x) + y = 1.2, 2x + cos(y) = 2")
		fmt.Println("2. x^2 + y^2 = 4, exp(x) + y = 1")
		fmt.Println("3. x + y^2 = 4, x^2 + y = 3")
		fmt.Println("4. sin(x+y) = 1.2x, x^2 + y^2 = 1")
		fmt.Println("5. x^2 - y = 1, x + y^2 = 1")
		fmt.Scan(&sysChoice)

		var systemFunc func(float64, float64) (float64, float64)
		var defaultX0, defaultY0 float64

		switch sysChoice {
		case 1:
			systemFunc = system1
			defaultX0, defaultY0 = 0.5, 0.5
		case 2:
			systemFunc = system2
			defaultX0, defaultY0 = 1.0, 1.0
		case 3:
			systemFunc = system3
			defaultX0, defaultY0 = 1.0, 1.0
		case 4:
			systemFunc = system4
			defaultX0, defaultY0 = 0.5, 0.5
		case 5:
			systemFunc = system5
			defaultX0, defaultY0 = 1.0, 1.0
		default:
			fmt.Println("Неверный выбор системы.")
			return
		}

		err := plotSystem(
			func(x, y float64) float64 { f, _ := systemFunc(x, y); return f },
			func(x, y float64) float64 { _, f := systemFunc(x, y); return f },
			defaultX0-4, defaultX0+4,
			defaultY0-4, defaultY0+4,
			0, 0,
			"system_preview.png",
		)
		if err != nil {
			fmt.Println("Ошибка при построении графика системы:", err)
		} else {
			fmt.Println("График системы сохранен в файл system_preview.png")
			fmt.Println("Пожалуйста, посмотрите график для выбора начальных приближений")
		}

		var x0, y0, epsilon float64
		fmt.Printf("Введите начальные приближения x0 и y0 (рекомендуемые около %v, %v):\n", defaultX0, defaultY0)
		fmt.Scan(&x0, &y0)

		fmt.Println("Введите погрешность вычисления (epsilon):")
		fmt.Scan(&epsilon)

		phi1 := func(x, y float64) float64 {
			f1, _ := systemFunc(x, y)
			return x - 0.1*f1
		}
		phi2 := func(x, y float64) float64 {
			_, f2 := systemFunc(x, y)
			return y - 0.1*f2
		}

		x, y, iterations, errors, converged := simpleIterationSystem(phi1, phi2, x0, y0, epsilon, 1000)
		if converged {
			fmt.Printf("\nРешение системы: x = %v, y = %v\n", x, y)
			fmt.Printf("Количество итераций: %d\n", iterations)
			fmt.Printf("Вектор погрешностей: %v\n", errors)
		} else {
			fmt.Println("Метод не сошелся за указанное число итераций")
		}

		f1, f2 := systemFunc(x, y)
		fmt.Printf("Проверка решения: f1(x, y) = %v, f2(x, y) = %v\n", f1, f2)

		file, err := os.Create("results_system.txt")
		if err != nil {
			fmt.Println("Ошибка при создании файла:", err)
			return
		}
		defer file.Close()

		fmt.Fprintf(file, "Решение системы: x = %v, y = %v\n", x, y)
		fmt.Fprintf(file, "Количество итераций: %d\n", iterations)
		fmt.Fprintf(file, "Вектор погрешностей: %v\n", errors)
		fmt.Fprintf(file, "Проверка решения: f1(x, y) = %v, f2(x, y) = %v\n", f1, f2)
	} else {
		fmt.Println("Неверный выбор задачи.")
	}
}

// выбор x_0 и достаточное условие сходимости в ньютоне
