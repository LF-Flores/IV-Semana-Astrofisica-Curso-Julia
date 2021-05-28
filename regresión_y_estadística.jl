### A Pluto.jl notebook ###
# v0.14.7

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

# ╔═╡ f032de0c-a3c4-48f8-9b95-4227c490a022
begin
	using Pkg; Pkg.activate(".")
	using PlutoUI, Distributions, StatsPlots, KernelDensity
end

# ╔═╡ f516adbd-33f2-4b89-a0e9-daa8e2e09bc2
PlutoUI.TableOfContents(title = "Contenido 💻")

# ╔═╡ ec5126de-bc13-11eb-39a4-f136d4015cf4
md"
# Regresión y estadística paramétrica
## Repaso de estadística básica

La estadística es una disciplina aplicada que le interesa la colección, organización, análisis, interpretación y presentación de los datos. Como principal auxiliar matemático utiliza a la Probabilidad, un subcampo de la teoría de la medida (que a su vez podría considerarse subcampo del análisis). 

Comenzaremos repasando probabilidad básica. Una buena intuición inicial de la probabilidad viene una cita de Henri Poincaré: *'El azar no es más que la medida de nuestra ignorancia'*.

Pero la ignorancia de un fenómeno puede significar muchas cosas. No es lo mismo aceptar ignorancia de cómo las partículas de un gas se mueven a aceptar ignorancia de cuáles son sus velocidades individuales.

Cada aspecto que ignoramos de un fenómeno debe ser reemplazado por el tipo adecuado de 'azar'. A este reemplazo adecuado le llamamos variable aleatoria: Una función que imaginamos que para cualquier evento que pueda subyacer al fenómeno, nos asigna *algún valor resumen* tal que no tendremos necesariamente claro cuál valor le corresponde a cada configuración del sistema pero tendremos información sobre cuáles valores son más probables que otros.

Como ejemplo, sea $X$ una variable aleatoria que represente la ignorancia aceptada sobre la dinámica de un dado siendo lanzado. Independiente de sus ecuaciones del movimiento, el estado final de un dado siendo lanzado puede resumirse en la cara que ha sido expuesta hacia arriba una vez que yazca en el suelo.

Por ello, decimos que tenemos 6 **consecuencias** de nuestros eventos (imágenes de $X$ sobre todo $\Omega$), $\{1, 2, 3, 4, 5, 6\}$, donde **evento** se debe entender como cualquier configuración en el espacio de estados del sistema dinámico, el cual admitimos ignorar en detalle (es decir, las preimágenes de $X$ que se resumen por ella). Es decir, el dominio de $X$ es $\Omega$ mientras que su rango es $\{1, 2, 3, 4, 5, 6\}$, siendo $\mathbb{R}$ siempre su codominio.

A éste conjunto de eventos le llamamos *espacio muestral*. Como comentado por Poincaré, el azar es una **medida de nuestra ignorancia** y la forma en que lo medimos es mediante la probabilidad. 

Sabemos que si lanzamos el dado, *alguno* de los 6 estados debió haber ocurrido. En ese sentido, decimos que el conjunto entero $\Omega$ tiene medida del 100% (acorde a las consecuencias que podemos tener) o bien, 1 denotando la enteridad. Cualquier subconjunto de $\Omega$ entonces debe tener medida menor o igual a 1.

Sabemos además que, en este caso, al menos un resultado debe ocurrir, por lo que el conjunto $\{n\}$ para cualquier $n \in \{1,2,3,4,5,6\}$ debe poder ser producido por al menos algún evento de $\Omega$. Consideremos el subconjunto $\{\omega \in \Omega | X(\omega) = n\}$. Este subconjunto se suele resumir por $[X = \omega]$ y es el conjunto de todos los eventos (aunque no podamos expresarlos explícitamente) que resultan en el valor resumen $n$ tras ser 'visto mediante el lente de $X$'.

Lo anterior puede no ser siempre cierto, pues existen sistemas dinámicos para los cuales una variable aleatoria que resuma su dinámica resulte tener eventos que nunca suceden. La noción de subconjuntos de $\Omega$ (eventos) que tienen medida $0$ es un tecnicismo que surge de la necesidad de definir apropiadamente *cuáles* subconjuntos de $\Omega$ deberían ser *medibles* y cuáles no, a pesar de que su medida sea $0$. El conjunto de los subconjuntos de $\Omega$ que **sí son medibles** se llama $\sigma$-álgebra y es usualmente denotada por $\Sigma$.

En todo caso, toda la información de cómo aceptamos nuestra ignorancia sobre este fenómeno está en una función $P: \Sigma \rightarrow [0,1]$ que llamamos **función de medida de probabilidad**. Ésta asigna a cada evento una medida que es interpretada como la probabilidad (o frecuencia relativa al resto de los eventos) con la que ocurre.

Recordando que la variable aleatoria es una función de la forma: $X: \Omega \rightarrow \mathbb{R}$, resulta que podemos siempre expresar la enteridad de $\Omega$ mediante la siguiente familia de conjuntos: 

$\{\omega \in \Omega | X(\omega) \leq x\} \quad \forall x \in \mathbb{R}$

Y es conveniente definir una función, llamada función de probabilidad acumulada, como:

$F_X(x) := P(\{\omega \in \Omega | X(\omega) \leq x\}) = P([X \leq x]) = P(X \leq x)$.

En el área de estadística, la mayoría de distribuciones estudiadas tienen la propiedad adicional que esta función $F_X: \mathbb{R} \rightarrow \mathbb{R}$ es derivable y definimos a:

$f_X(x) = \frac{d}{dx} F_X(x)$

la cual es llamada **función de densidad de probabilidad** o solamente **densidad** y a veces denotada solo por $f(x)$ cuando se sobre entiende la variable aleatoria en cuestión.

observemos algunas *distribuciones de probabilidad* (refiriéndonos a variables aleatorias que siguen una función de distribución específica) populares:
"

# ╔═╡ d4088ab5-c0e0-4b48-b778-f2ef1ee96c06
md"
## Distribuciones paramétricas populares

### La distribución uniforme

En Julia tenemos la distribución uniforme por defecto cuando se corre la función `rand`:"

# ╔═╡ 4c6b6aeb-95ec-436f-b4fb-2d02ebffcb76
rand() # Esto es lo que llamamos la distribución U([0,1])

# ╔═╡ eacc76f5-725f-4932-8e7c-500359295793
rand(10)

# ╔═╡ 66861f13-e67b-4065-9950-b49f202a08aa
md"números de tiros = $(@bind número_de_tiros Slider(0:50:2000, show_value = true))"

# ╔═╡ 39e512bd-e87e-47a0-93b0-277342c49a56
tiros_de_dado = rand(1:6, número_de_tiros);

# ╔═╡ 4db884d4-800d-4c9b-8c81-42e7b667ccfd
histogram(tiros_de_dado, bins = range(1, 6, step = 1))

# ╔═╡ 3b79fe4f-71ba-410d-adae-7ec46a599032
md" ### La distribución Normal"

# ╔═╡ 30758af4-dd46-4aa0-b7d8-93d368591396
d₁ = Normal()

# ╔═╡ 3b79dcc3-c7b9-4108-b35a-195bf1c962c4
md"n = $(@bind n₁ Slider(20:100:90000, show_value = true))"

# ╔═╡ 9a6a8750-ce42-41ad-98a4-e2b8f707bb5b
sample₁ = rand(d₁, n₁);

# ╔═╡ 31b7baaa-1f86-46e8-b283-6b1f2ee40b06
fit₁ = fit(Normal, sample₁)

# ╔═╡ 3e078961-eac6-479b-bee4-d361c74e4804
fit₁.μ, fit₁.σ - 1 

# ╔═╡ fe617c37-bd4e-4bf9-932f-cbaac349f187
mean(sample₁), std(sample₁) - 1

# ╔═╡ f3eca08a-3f59-4ab9-bf7a-3df35294032f
density(sample₁)

# ╔═╡ bc35e4ad-98b3-490f-a405-ece1a86c76c0
md" ### La distribución Chi"

# ╔═╡ 553719e5-2dde-40ba-85a1-56f4bf1a55bc
d₂ = Chi(3)

# ╔═╡ c24fd5d4-8abe-43f0-bb61-232264185d58
μᵪ = 2sqrt(2)/sqrt(π)

# ╔═╡ da85cab1-8721-4dcd-bc33-117da8124c82
md"n = $(@bind n₂ Slider(20:100:90000, show_value = true))"

# ╔═╡ 49a82540-6be5-4779-89be-c587bab7e133
sample₂ = rand(d₂, n₂);

# ╔═╡ eec2243b-6675-4836-9ea9-1dd38daf0511
mean(sample₂) - μᵪ, std(sample₂) - sqrt((3 - μᵪ^2))

# ╔═╡ d8e2f0cc-4721-47d6-aeb0-1bfb4bdf0499
density(sample₂)

# ╔═╡ 6d5afab8-2aca-4445-9d6c-f935a94f4f81
md"
## La distribución Gamma
"

# ╔═╡ 0ae3a70a-3ea5-448e-9ebb-8083e8067b3f
md"T = $(@bind T Slider(1:50:500, show_value = true))"

# ╔═╡ 38110dcc-e84d-4447-96bd-278b876f47be
d₃ = Gamma(3/2, T)

# ╔═╡ 8fea8b30-945d-4c24-a7d7-48c7f568a9af
sample₃ = rand(d₃, 10000);

# ╔═╡ e8491c29-163f-4230-9848-c205809562ee
density(sample₃)

# ╔═╡ 33bd3752-ece1-4acb-8a2c-c34f3fa85207
md"## Más distribuciones

Podemos encontrar más distribuciones en la [**documentación de** `Distributions.jl`](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Gamma), donde muestran todas las distribuciones continuas, discretas y truncadas que manejan, tanto a nivel univariado como multivariado y que sirven para modelar diversidad de fenómenos físicos para los cuales aceptamos un grado de estocástidad debido a ignorancia que se desea aceptar en los detalles de su dinámica.

Por una parte, podemos ver un ejemplo más práctico, en donde podemos hablar de modelos mixtos, los cuales surgen cuando nuestra data es resultado de varios procesos/dinámicas cuyos efectos son igualmente dominantes."

# ╔═╡ ef5d73c0-4b8c-4185-8e6a-d386dfa9d808
test = MixtureModel(Normal[
   Normal(-2.0, 1.2),
   Normal(0.0, 1.0),
   Normal(3.0, 2.5)])

# ╔═╡ b130f24b-a421-433b-be5e-8496a0117cae
md"c₁ = $(@bind c₁ Slider(1:0.1:10, show_value = true))

c₂ = $(@bind c₂ Slider(1:0.1:10, show_value = true))

c₃ = $(@bind c₃ Slider(1:0.1:10, show_value = true))"

# ╔═╡ 4fc7016d-53de-4d4c-8151-e8b424892efd
density(rand(MixtureModel(map(u -> Normal(u, 1.5), [-5, 0, 4]), [c₁, c₂, c₃]./(c₁ + c₂ + c₃)), 10000))

# ╔═╡ cf5ddf0f-8136-43d2-89d3-b422736cd553
md"k₁ = $(@bind k₁ Slider(1:0.1:10, show_value = true))

k₂ = $(@bind k₂ Slider(1:0.1:10, show_value = true))
"

# ╔═╡ dcd2ae59-872f-4dca-8ead-5aff5dd6a961
density(rand(MixtureModel([Exponential(1), Normal(2,0.5)], [k₁, k₂]./(k₁ + k₂)), 10000))

# ╔═╡ 3036f2af-d5f8-482f-82d2-fbbd9ce61c2b
md"## Estimación del kernel de una densidad de probabilidad"

# ╔═╡ 86476b19-04a4-40d4-b71f-af3b2d814012
data = rand(MixtureModel(map(u -> Normal(u, 1.5), [-5, 0, 4]), [0.496552, 0.358621, 1 - 0.496552 - 0.358621]), 10000) 

# ╔═╡ 8acb86ac-aaaa-489c-b95b-edc67bc75507
kernel_estimado = kde(data)

# ╔═╡ 7753565a-7d21-49fe-ad77-052b1e75997a
begin
	density(data, label = "Densidad original")
	plot!(kernel_estimado, linestyle = :dash, linewidth = 4, label = "Densidad Estimada")
	title!("Estimación con KDE con filtro Gaussiano")
end

# ╔═╡ Cell order:
# ╠═f032de0c-a3c4-48f8-9b95-4227c490a022
# ╟─f516adbd-33f2-4b89-a0e9-daa8e2e09bc2
# ╟─ec5126de-bc13-11eb-39a4-f136d4015cf4
# ╟─d4088ab5-c0e0-4b48-b778-f2ef1ee96c06
# ╠═4c6b6aeb-95ec-436f-b4fb-2d02ebffcb76
# ╠═eacc76f5-725f-4932-8e7c-500359295793
# ╟─66861f13-e67b-4065-9950-b49f202a08aa
# ╠═39e512bd-e87e-47a0-93b0-277342c49a56
# ╠═4db884d4-800d-4c9b-8c81-42e7b667ccfd
# ╟─3b79fe4f-71ba-410d-adae-7ec46a599032
# ╠═30758af4-dd46-4aa0-b7d8-93d368591396
# ╠═9a6a8750-ce42-41ad-98a4-e2b8f707bb5b
# ╠═31b7baaa-1f86-46e8-b283-6b1f2ee40b06
# ╠═3e078961-eac6-479b-bee4-d361c74e4804
# ╠═fe617c37-bd4e-4bf9-932f-cbaac349f187
# ╟─3b79dcc3-c7b9-4108-b35a-195bf1c962c4
# ╠═f3eca08a-3f59-4ab9-bf7a-3df35294032f
# ╟─bc35e4ad-98b3-490f-a405-ece1a86c76c0
# ╠═553719e5-2dde-40ba-85a1-56f4bf1a55bc
# ╠═49a82540-6be5-4779-89be-c587bab7e133
# ╠═c24fd5d4-8abe-43f0-bb61-232264185d58
# ╠═eec2243b-6675-4836-9ea9-1dd38daf0511
# ╟─da85cab1-8721-4dcd-bc33-117da8124c82
# ╠═d8e2f0cc-4721-47d6-aeb0-1bfb4bdf0499
# ╟─6d5afab8-2aca-4445-9d6c-f935a94f4f81
# ╠═38110dcc-e84d-4447-96bd-278b876f47be
# ╠═8fea8b30-945d-4c24-a7d7-48c7f568a9af
# ╟─0ae3a70a-3ea5-448e-9ebb-8083e8067b3f
# ╠═e8491c29-163f-4230-9848-c205809562ee
# ╟─33bd3752-ece1-4acb-8a2c-c34f3fa85207
# ╠═ef5d73c0-4b8c-4185-8e6a-d386dfa9d808
# ╟─b130f24b-a421-433b-be5e-8496a0117cae
# ╠═4fc7016d-53de-4d4c-8151-e8b424892efd
# ╟─cf5ddf0f-8136-43d2-89d3-b422736cd553
# ╠═dcd2ae59-872f-4dca-8ead-5aff5dd6a961
# ╟─3036f2af-d5f8-482f-82d2-fbbd9ce61c2b
# ╠═86476b19-04a4-40d4-b71f-af3b2d814012
# ╠═8acb86ac-aaaa-489c-b95b-edc67bc75507
# ╠═7753565a-7d21-49fe-ad77-052b1e75997a
