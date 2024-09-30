# stan-intro-gab24

## Deconstruir el modelado no eStan complicado. Una introducción a Stan.

Taller para la reunión del Grupo Argentino de Bioestadística 2024
(Buenos Aires).

Descripción:

Si la variable explicatoria es categórica, ANOVA; si es cuantitativa,
regresión; si no se cumple el supuesto de normalidad, GLM, y así hasta
las ramas más recónditas de la clave dicotómica. A veces la estadística
parece una bolsa de métodos inconexos acompañada de una receta para
usarlos. Bajo este enfoque, elegimos qué modelo utilizar siguiendo la
receta y evaluando la disponibilidad de software. Pero las herramientas
computacionales actuales permiten adoptar otro enfoque: partiendo de
nuestras preguntas u objetivos, podemos formular un modelo que
represente el proceso bajo estudio de manera razonable, y luego buscar
la forma de estimarlo. En este taller veremos cómo dar los primeros
pasos hacia un modelado más libre de recetas usando Stan, una plataforma
estadística orientada a la estimación de modelos estadísticos con
enfoque Bayesiano.

Usaremos Stan desde R, en donde también estimaremos modelos maximizando
la likelihood. Para seguir el taller se requieren los siguientes
paquetes instalados: 
- rstan
- bbmle
- DHARMa
- tidyverse
- ggdensity
- viridis.

Se recomienda seguir el taller usando el archivo `modelos.Rmd`. En caso de usar
RStudio, lo ideal es tener una versión $\ge$ 1.4 para usar el modo de
edición "visual", así el texto se diferencia mejor del código y las
ecuaciones se ven como se debe. En caso de no contar con los paquetes al día para ir corriendo el código, se puede abrir el archivo `modelos.html` en un navegador para ver el código desde cerca.
