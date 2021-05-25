### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ e9f41c8b-9400-4d2d-a2e0-db016b62b69c
begin
	using Pkg; Pkg.activate(".")
	using Plots, PlutoUI
end

# ╔═╡ ae6ae8e9-8d53-4a5a-929f-e02fa3a10901
PlutoUI.TableOfContents(title = "Contenidos 💻")

# ╔═╡ 3d577cb4-9418-45ee-a968-62f3b5287759
md"# Segundo día: Introducción al pensamiento computacional y algorítmico
## Minimización en una variable

considere la siguiente función:"

# ╔═╡ 8531e010-baf8-11eb-2b37-6be3738f5f8e
f(x) = (x^2 - 1)*(x-2)

# ╔═╡ cef9b84e-00e6-4c8f-9694-03b2a8557eb4
plot(f, -1.5:0.1:2.7) 

# ╔═╡ ee1fe570-b26e-4adc-a5f7-4fbc77c5a29f
md"Esta función, al ser un polinomio, es fácil encontrar sus mínimos locales y mínimo global utilizando el cálculo. No obstante, podemos usarlo de modelo de alguna función con características similares pero forma analítica mucho más compleja, o potencialmente una función empírica para la cual solamente tengamos los datos obtenidos de observaciones.

### Encontrando un intervalo inicial

El primer paso para encontrar el mínimo local de una función general unimodal (que tiene solo un mínimo) es encontrar un intervalo $[a,b]$ tal que el punto minimizante $x_{\text{mín}}$ esté en $(a,b)$. En caso de funciones que no son unimodales, podemos tener el riesgo de encontrar un mínimo local que no sea un mínimo global. 

Ahora, incluso en funciones unimodales, debemos encontrar un intervalo $[a,b]$ válido como descrito anteriormente pues puede ser el caso que comencemos la búsqueda del mínimo erroneamente con un intervalo $[a,b]$ tal que $x < a < b$ o $a < b < x$ para los $x$ verificados.

Un método para poder encontrar un bracket inicial con el cual trabajar es el siguiente:

1. Comenzamos con un punto inicial cualquiera, $x$. Entre más cerca de $x_{\text{mín}}$ estemos, mejor (aunque a veces no hay forma de saber).
2. Dicho punto y su imagen: $x, f(x)$, van a ser respectivamente la primera propuesta de $a$ (en $[a,b]$) y su imagen, $f(a)$.
3. Tomamos un pequeño paso a la izquierda (digamos de tamaño $s$) para obtener la primera propuesta de $b$ y su imagen, $f(b)$.
4. Puede ser el caso que nuestro paso de tamaño $s$ nos haya llevado a un punto más alto y no uno más bajo (¡recordemos que deseamos minimizar la función!). En ese caso, debemos intercambiar los roles de $a$ con $b$ y los de sus imágenes. 
5. Una vez que sepamos que estamos en dirección de descenso, debemos seguir tomando pasos (si deseamos, siempre de tamaño $s$ o, en general, de un nuevo tamaño $k \times s$ cada vez). 
6. El algoritmo de detendrá en el momento que la imagen del punto actualizado, $x_{\text{prev}} + k\times s$, sea mayor a la imagen de $x_{\text{prev}}$ ya que, bajo la premisa de que ya estábamos en dirección de descenso, esto solo puede significar que hemos saltado sobre el punto mínimo y comenzado a subir. 
7. En conclusión de ello, los valores actuales en esa iteración para $a$ y $b$ serán los que garanticen que $x_{\text{mín}} \in [a,b]$.

Implementemos el algoritmo a continuación:
"

# ╔═╡ 69cdb427-8822-4cd6-9be7-ee4593ad7c21
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
	a, ya = x, f(x)
	b, yb = a + s, f(a + s)
	if yb > ya
		a, b = b, a
		ya, yb = yb, ya
		s = -s
	end
	while true
		c, yc = b + s, f(b + s)
		if yc > yb
			return a < c ? (a, c) : (c, a)
		end
		a, ya, b, yb = b, yb, c, yc
		s *= k
	end
end		

# ╔═╡ 61db35aa-8511-4241-8824-a0f3d9794ca8
md"Los parámetros $s$ y $k$ del algoritmo son llamados *hiperparámetros* ya que, más que parámetros que nos definen un comportamiento específico de nuestro algoritmo, son parámetros que tienen que ser decididos en torno al tipo de función, $f$, que estemos tratando, y tendrán mejor o peor rendimiento dependiendo de las cualidades de $f$.

Podemos ilustrar como la elección del intervalo con el mínimo depende de éstos dos para nuestra función actual:"

# ╔═╡ f92fdc19-1d54-492e-9690-f8a0bcadb425
md"s = $(@bind s Slider(0.02:0.02:1, show_value = true, default = 0.01))

   k = $(@bind k Slider(1.1:0.1:4.0, show_value = true, default = 2.0))"

# ╔═╡ 976e8f79-5b81-4ad9-af58-e84942ed7db5
let (a, b) = bracket_minimum(f, 1; s=s, k=k)
	plot(f, -1.5:0.1:2.7, label = "f(x)")
	plot!([(a, 0), (b, 0)], seriestype = :scatter, 
		  label = "intervalo encontrado", markershape = :vline,
	      markersize = 6)
end

# ╔═╡ 02f5dd5e-c221-4e21-98d4-958e8d7ac139
md"Una vez encontrado este intervalo, la siguiente tarea es encontrar el punto mínimo de la función, $x_{\text{mín}}$. Esto puede ser logrado por diversos algoritmos, pero en este caso deduciremos uno llamado **búsqueda de ajuste cuadrático** (o quadratic fit search).

### Búsqueda de ajuste cuadrático

Este deriva su nombre de su estrategia: Ajustar a $f$ en $[a,b]$ a una función cuadrática, bajo la premisa que muchas funciones se comportan similar a polinomios cuadráticos en vecindad de un mínimo. Esto quiere decir que, entre más pequeño y cercano sea $[a,b]$ del punto mínimo, mejor rendimiento tendremos.

Comencemos con la deducción:
1. Para fines de notación, llamemos a lo que previamente era el intervalo $[a,b]$ a $[a,c]$ ya que vamos a querer llevar rastro de un punto interior, que llamaremos $b$, tal que: $a < b < c$ para servir como propuesta de $x_{\text{mín}}$, la cual iremos actualizando.
2. Definamos nuestra función de ajuste cuadrático: $q(x) = p_1 + p_2x + p_3x^2$
3. Podemos hacer que $q(x)$ defina una cuadrática que pase por todos nuestros puntos conocidos: $(a, f(a)), (b, f(b)), (c, f(c))$ al encontrar los coeficientes $p_i$ adecuados. Como notación los $f(k)$ ahora será denotado $y_k$ para $k = a,b,c$. Expresamos las relaciones resultantes:

   - $y_a = p_1 + p_2a + p_3a^2$
   - $y_b = p_1 + p_2b + p_3b^2$
   - $y_c = p_1 + p_2c + p_3c^2$

Esto se puede resumir mejor en la siguiente ecuación matricial:

$\begin{bmatrix}
y_a \\ y_b \\ y_c
\end{bmatrix}
= \begin{bmatrix}
1 & a & a^2 \\ 1 & b & b^2  \\ 1 & c & c^2 
\end{bmatrix} 
\begin{bmatrix}
p_1 \\ p_2 \\ p_3 
\end{bmatrix}$

Para resolver esta ecuación basta con encontrar la matriz inversa, esto ya que podemos calcular que el determinante de esta matriz es: $-(a-b)(a-c)(b-c)$ (¡Verifíquenlo!) y ésta expresión solamente puede ser cero si al menos alguno de los tres puntos es el mismo. 

Asumiendo entonces que la matriz es no singular, tenemos:

$\begin{bmatrix}
p_1 \\ p_2 \\ p_3 
\end{bmatrix}
= \begin{bmatrix}
1 & a & a^2 \\ 1 & b & b^2  \\ 1 & c & c^2 
\end{bmatrix}^{-1} 
\begin{bmatrix}
y_a \\ y_b \\ y_c
\end{bmatrix}$

la cual puede ser resuelta con su método favorito. Ya que tenemos el determinante, podemos utlizar la [**regla de Cramer**](https://en.wikipedia.org/wiki/Cramer%27s_rule) para obtener que nuestro polinomio cuadrático ajustado es:

$q(x) = y_a\frac{(x-b)(x-c)}{(a-b)(a-c)} + y_b\frac{(x-a)(x-c)}{(b-a)(b-c)} + y_c\frac{(x-b)(x-a)}{(c-b)(c-a)}$

Calculando $q'(x)$, la derivada de $q(x)$, e igualando a cero podemos obtener el mínimo global de la cuadrática, el cual sirve como una buena propuesta para el punto mínimo $x_{\text{mín}}$ de nuestra función. Estos comenzarán a coincidir a medida el intervalo sea más pequeño y $f$ parezca cada vez más una cuadrática en el intervalo $[a,c]$. 

Motivados por ello, debemos proponer una forma de actualizar el intervalo $[a,c]$ utilizando el valor de la nueva propuesta óptima:

$x^* = \frac{1}{2} \frac{y_a(b^2-c^2) + y_b(c^2-a^2) + y_c(a^2-b^2)}{y_a(b-c) + y_b(c-a) + y_c(a-b)}$

que vendrá a substituir nuestra propuesta inicial $b$. Bosquejemos el algoritmo a continuación:

1. Definamos un número $n$ de iteraciones máximas que deseamos hacer. 

   Este $n$ podría estar controlado por algo más sofisticado como reconocer que nuestro algoritmo generará puntos en una sucesión convergente la cual también es de Cauchy y por ende podemos llevar rastro de la diferencia entre puntos propuestos consecutivos y saber que el resultado se irá haciendo más pequeno. Entonces, la condición de paro puede ser una tolerancia $\epsilon$ bajo la cual nos baste que yazca la diferencia de nuestras propuestas.


2. De nuestros tres puntos iniciales $a,b,c$ ($a$ y $c$ siendo los límites del intervalo que contiene el mínimo y $b$ un punto interior arbitrario), podemos calcular $x^*$ en base a la fórmula deducida anteriormente.


3. El $x^*$ encontrado puede estar a la izquierda o a la derecha de $b$ (pues, como será demostrable luego, como nuestro algoritmo converge, si $x^*$ coincide con $b$, quiere decir que $b$ ya era nuestro punto óptimo desde el inicio y podemos retornarlo).

   Independiente de dónde caiga, tendremos una estrategia general de actualización. Consideremos sin pérdida de generalidad que cae adelante de $b$. Entonces: $x^* \in (b, c]$ (en el caso que haya caído a la derecha, $x^* \in [a, b)$).


4. Verificamos si la imagen, $f(x^*)$, de nuestra propuesta es mayor o menor que $y_b$:
   - Si es mayor, quiere decir que nos hemos saltado a $x_{\text{mín}}$ y podemos utilizar a $x^*$ como nuestro nuevo $c$ para achicar el intervalo.
   - Si es menor, quiere decir que $x_{\text{mín}}$ aún no ha sido alcanzado (está por enfrente tanto de $x^*$ como de $b$) y entonces podemos descartar a $a$ como intervalo interior en favor de $b$.


5. Realizada esta actualización del intervalo, podremos repetir el paso 2 para encontrar una nueva propuesta y repetir los pasos 3 y 4 hasta cumplir nuestra condición de paro.


6. Una vez detenido el algoritmo, podemos devolver los últimos valores de los tres números: $a,b,c$ donde $b$ será nuestra mejor propuesta para $x_{\text{mín}}$ mientras que $[a,c]$ seguirá siendo un intervalo que garantice contener a $x_{\text{mín}}$.

Prosigamos a implementarlo:
"

# ╔═╡ bc48831c-1c0b-4091-b8b5-07384de3039b
function quadratic_fit_search(f, a, b, c, n)
	ya, yb, yc = f(a), f(b), f(c)
	for i in 1:n
		x = 0.5*(ya*(b^2-c^2)+yb*(c^2-a^2)+yc*(a^2-b^2)) /
		(ya*(b-c)
		+yb*(c-a)
		+yc*(a-b))
		yx = f(x)
		if x > b
			if yx > yb
				c, yc = x, yx
			else
				a, ya, b, yb = b, yb, x, yx
			end
		elseif x < b
			if yx > yb
				a, ya = x, yx
			else
				c, yc, b, yb = b, yb, x, yx
			end
		end
	end
	return (a, b, c)
end

# ╔═╡ ad161424-fcd7-4d3e-8d65-c1d9be8ab92f
md"Quisieramos poder juntar ambos algoritmos: El de la búsqueda del intervalo con mínimo y la búsqueda de ajuste cuadrático. Para ello, aprovechamos el **multiple dispatch** de Julia para definir otro método de la función `quadratic_fit_search` que esta vez solo tome a $f$, a un punto inicial $x$ para la búsqueda del intervalo y el número de iteraciones $n$.

Para ello, decidimos que nuestra propuesta inicial, $b \in [a,b]$ siempre será el punto medio del intervalo. Prosegimos a la implementación: "

# ╔═╡ 701dbc97-ad20-4a22-9dfc-8bfc2e4b806f
function quadratic_fit_search(f, x, n)
	a, c = bracket_minimum(f, x)
	b = (a+c)/2
	quadratic_fit_search(f, a, b, c, n)
end

# ╔═╡ 87292b27-8aef-495c-86f9-fca65993c382
md"Exploremos ahora cómo la optimización de nuestra función $f$ depende de la elección de $x$ inicial y del número de iteraciones."

# ╔═╡ b4cd0cc3-bafe-45ab-b549-8a1d1e471195
óptimo_real = (2+sqrt(7))/3

# ╔═╡ 3475b4a0-d3b9-48c3-8cd5-70afad0e0aae
md"x = $(@bind x Slider(0.5:0.1:3, show_value = true))

   n = $(@bind n Slider(1:1:10, show_value = true))"

# ╔═╡ 744e4de6-d4fe-4778-9291-ac8b66c887e4
bracket_minimum(f, x)

# ╔═╡ 269f3e2f-df18-42a6-85de-1bfa9f84b358
a_opt, x_opt, c_opt = quadratic_fit_search(f, x, n)

# ╔═╡ f1fcf108-bd8f-44d8-b28c-00228db3d092
x_opt

# ╔═╡ 71722f4d-e8a9-4758-ab1d-477e5bf730c3
óptimo_real - x_opt

# ╔═╡ 87160279-33c2-40ca-908b-37dad0785d80
mean_error = sum(óptimo_real .- (x -> quadratic_fit_search(f, x, n)[2]).(0.5:0.1:3))/length(0.5:0.1:3)

# ╔═╡ bd2dca3a-dce7-4431-a5b3-a736b89d5493
begin
	plot(f, 0.5:0.01:2.2, label = "f(x)") 
	plot!([(óptimo_real, f(óptimo_real))], seriestype = :scatter, label = "óptimo real")
	plot!([(x_opt, f(x_opt))], seriestype = :scatter, label = "punto predicho")
	plot!([(a_opt, -0.65), (c_opt, -0.65)], seriestype = :scatter, 
		  label = "intervalo final", markershape = :vline,
	      markersize = 10)
end

# ╔═╡ ec124952-0a0d-4d24-b0b7-556945dea6f7
md"
## Minimización multivariable

El tema de minimización en una variable es muchísimo más extenso de lo que discutimos. Esto significa por supuesto que veremos muy por encima la minimización multivariable, pero el propósito de presentar al menos un método es que sea más evidente la importancia de construir sobre fundamentos rigurosos hacia métodos que son utilizados en algoritmos muy importantes de la industria y frontera de la ciencia.

El algoritmo que visitaremos se llama **line search** y es parte de una familia de algoritmos utilizados para encontrar una propuesta de un punto más cercano al mínimo de alguna función $f: \mathbb{R}^n \rightarrow \mathbb{R}$ tras ya conocer una *dirección óptima de descenso* (la cual no discutiremos cómo encontrar, pues involucra otra familia de métodos llamados descenso del gradiente, la cual no podremos discutir por tiempo).

Consideremos la siguiente función:" 

# ╔═╡ 1e067699-f725-4f45-9ba4-fafb315ee8cd
g(x, y) = 2*x^2 + y^2

# ╔═╡ f9434958-f02d-484d-a54c-1039f5f91949
begin
	xs = as = collect(-1:0.1:1.0)
	x_grid = [x for x = xs for y = as]
	a_grid = [y for x = xs for y = as]
	plot(x_grid, a_grid, g.(x_grid, a_grid), st = :surface, xlabel = "longer xlabel", ylabel = "longer ylabel", zlabel = "longer zlabel")
end

# ╔═╡ 1a10f9e8-044d-4b48-8887-854a029f25d4
md"
Imaginemos que conocemos un vector $\mathbf{d} \in \mathbb{R}^2$ que define una dirección de desenco desde un punto inicial, $\mathbf{x_0}\in \mathbb{R}^2$. Esta dirección de descenso es necesaria de saber similar a cómo anteriormente ocupamos verificar en `bracket_minimum` que nuestros pasos `s` no nos llevaran a un valor mayor y por ende más lejos del mínimo.

Este valor de $\mathbf{d}$ es posible e intuitivo de definir como la dirección normalizada opuesta al gradiente (ya que el gradiente es la dirección del cambio más inclinado positivo). Por ahora, lo dejaremos fijo a un valor y nos preocuparemos por bosquejar el algoritmo básico de *line search*:

1. Conociendo un punto inicial $\mathbf{x_0}$ y la siguiente dirección de descenso $\mathbf{d}$, podemos estudiar la **sección de la funcíón** $f$ definida por:

   $f_{\mathbf{d}} : \alpha \longmapsto f(\mathbf{x} + \alpha \mathbf{d})$

   Es decir, aumentar la preimagen de f a lo largo de la dirección $\mathbf{d}$. Como esta nueva función es de una dimensión, podemos encontrar el $\alpha_{\text{mín}}$ con los métodos de una variable.


2. Elegimos usar `quadratic_fit_search` junto con `bracket_minimum` para optimizar la función $f_\mathbf{d}$, aunque pudo haber sido cualquier otro método de optimización de una variable.


3. Retornamos como propuesta a un valor más cercano a $\mathbf{x}_\text{mín}$ al valor avanzado por el $\alpha$ \optimo: 

$\mathbf{x_1} = \mathbf{x_0} + \alpha_\text{mín} \mathbf{d}$

El algoritmo se espera iterar subsecuentes veces, con diferentes valores de $\mathbf{d}$ cada vez, y tomando como valores iniciales los \mathbf{x_k} que vayamos generando en las iteraciones anteriores hasta converger al mínimo de la $f$ multivariable.

Implementemos `line_search`:
"

# ╔═╡ 6058cd96-9a01-4e2e-a532-88ff48fab5ff
function line_search(f, x, d)
	objective = α -> f(x + α*d...)
	_, α, _ = quadratic_fit_search(objective, x'*d, 5)
	return x + α*d
end

# ╔═╡ 7ff4af40-f4d5-41c1-88e6-6add4870451d
md"Podemos verlo en acción en las siguientes celdas, notando que debemos aún colocar *'a mano'* el vector direccional:"

# ╔═╡ ed53175a-a94f-4149-8974-91c9f6699943
x₀ = [1,-1]

# ╔═╡ 558c4f2b-a484-4d3e-a1a1-f1cd4e6a035f
x₁ = line_search(g, x₀, [1/sqrt(2), 1/sqrt(2)])

# ╔═╡ b6c24326-fdd0-423b-8009-3dacad259ab5
x₂ = line_search(g, x₁, [1/sqrt(2), -1/sqrt(2)])

# ╔═╡ 7684c089-d99d-4221-ad9b-ae25c13f2708
x₃ = line_search(g, x₂, [1/sqrt(2), 1/sqrt(2)])

# ╔═╡ fec3f7de-07b0-401d-ba88-98aa47c398e4
line_search(g, x₃, [1/sqrt(2), -1/sqrt(2)])

# ╔═╡ c33963ba-f967-46b4-915b-22b4492f610e
md"
## Tarea 😠

**Intrucciones:** Debe hacer AL MENOS uno de los siguientes ejercicios, pero se recomienda realizarlos todos para poder practicar los conceptos vistos.

1. Verificar la deducción del punto $x^*$ en el algoritmo de la búsqueda de ajuste cuadrático con detalle.


2. Demostrar que en el algoritmo de la búsqueda de ajuste cuadrático, nuestra sucesión de propuestas: $\{x^*_k\}_{k=1}^\infty$ realmente converge al mínimo local $x_{\text{mín}}$. *Pista:* pueden utilizar el [teorema del sándwich](https://en.wikipedia.org/wiki/Squeeze_theorem) o probar que es una [sucesión de Cauchy](https://en.wikipedia.org/wiki/Cauchy_sequence).


3. Implementar una función `quadratic_fit_search(f, x, ϵ = 1e-3)` que continúe iterando sobre el quatratic search hasta que dos propuestas consecutivas tengan una diferencia menor a `ϵ`. note que `ϵ` es un parámetro opcional. 

   De esta manera podremos llamar `quadratic_fit_search(f, x)` para algún $x$ que creamos que está cerca del mínimo real y el algoritmo se encargue del resto sin preocuparnos de especificar todos los hiperparámetros. **Nota:** También es válido poner al resto de hiperparámetros como opcionales.

4. El `line_search` implementado es el más básico de toda la familia de métodos de line search. Diseñemos uno ligeramente más complicado al reconocer que el `line_search` normal puede ser costoso de ejecutar ya que en cada paso progresivo hacia el mínimo, debemos ejecutar una optimización completa de `quadratic_fit_search`. 

   Esto puede no ser necesario. Basta con que la función $f$ decrezca *lo suficiente*. ¿Cuánto es suficiente? Consideremos la siguiente condición:

   $f(\mathbf{x_{k+1}}) \leq f(\mathbf{x_{k}}) + \beta \alpha \nabla_{\mathbf{d_{k}}} f(\mathbf{x_{k}})$ 

   Un line search que se detiene cuando la condición anterior se cumple se conoce como *line search aproximado con condición de descenso suficiente* o *backtracking line search* y tiene como intuición subyacente que: $\nabla_{\mathbf{d_{k}}} f(\mathbf{x_{k}})$ es la pendiente inmediata de $f$ en el punto $\mathbf{x_{k}}$ y en dirección $\mathbf{d_{k}}$ y ésta define un descenso que haríamos a primer orden de aproximación, lo cual debería ser suficiente para detener la iteración. 

   No obstante, tenemos un hiperparámetro $\beta$ que controla cuánto porcentaje de ese descenso de primer orden consideramos suficiente. Por ejemplo, si $\beta = 0.5$, estaríamos considerando que la mitad de lo que el gradiente los dice que es suficiente descender es ya suficiente para nosotros y en caso que $\beta = 0$, cualquier descenso es suficiente.

   Implementa aquí el algoritmo de `backtracking_line_search(f, ∇f, x, d, α, β)` que con cualquier técnica de minimización de una dimensión logre minimizar una función unimodal multivariada más eficientemente que nuestra implementación de `line_search`.
"

# ╔═╡ Cell order:
# ╠═e9f41c8b-9400-4d2d-a2e0-db016b62b69c
# ╟─ae6ae8e9-8d53-4a5a-929f-e02fa3a10901
# ╟─3d577cb4-9418-45ee-a968-62f3b5287759
# ╠═8531e010-baf8-11eb-2b37-6be3738f5f8e
# ╠═cef9b84e-00e6-4c8f-9694-03b2a8557eb4
# ╟─ee1fe570-b26e-4adc-a5f7-4fbc77c5a29f
# ╠═69cdb427-8822-4cd6-9be7-ee4593ad7c21
# ╟─61db35aa-8511-4241-8824-a0f3d9794ca8
# ╟─f92fdc19-1d54-492e-9690-f8a0bcadb425
# ╟─976e8f79-5b81-4ad9-af58-e84942ed7db5
# ╟─02f5dd5e-c221-4e21-98d4-958e8d7ac139
# ╠═bc48831c-1c0b-4091-b8b5-07384de3039b
# ╟─ad161424-fcd7-4d3e-8d65-c1d9be8ab92f
# ╠═701dbc97-ad20-4a22-9dfc-8bfc2e4b806f
# ╟─87292b27-8aef-495c-86f9-fca65993c382
# ╠═744e4de6-d4fe-4778-9291-ac8b66c887e4
# ╠═269f3e2f-df18-42a6-85de-1bfa9f84b358
# ╠═f1fcf108-bd8f-44d8-b28c-00228db3d092
# ╠═b4cd0cc3-bafe-45ab-b549-8a1d1e471195
# ╠═71722f4d-e8a9-4758-ab1d-477e5bf730c3
# ╟─87160279-33c2-40ca-908b-37dad0785d80
# ╟─3475b4a0-d3b9-48c3-8cd5-70afad0e0aae
# ╟─bd2dca3a-dce7-4431-a5b3-a736b89d5493
# ╟─ec124952-0a0d-4d24-b0b7-556945dea6f7
# ╠═1e067699-f725-4f45-9ba4-fafb315ee8cd
# ╟─f9434958-f02d-484d-a54c-1039f5f91949
# ╟─1a10f9e8-044d-4b48-8887-854a029f25d4
# ╠═6058cd96-9a01-4e2e-a532-88ff48fab5ff
# ╟─7ff4af40-f4d5-41c1-88e6-6add4870451d
# ╠═ed53175a-a94f-4149-8974-91c9f6699943
# ╠═558c4f2b-a484-4d3e-a1a1-f1cd4e6a035f
# ╠═b6c24326-fdd0-423b-8009-3dacad259ab5
# ╠═7684c089-d99d-4221-ad9b-ae25c13f2708
# ╠═fec3f7de-07b0-401d-ba88-98aa47c398e4
# ╟─c33963ba-f967-46b4-915b-22b4492f610e
