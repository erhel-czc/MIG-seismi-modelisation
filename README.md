# MIG-seismi-modelisation

_Objectif : modéliser la sismicité induite par injection de fluide en géothermie en Alsace._

## Jour 1
- création d'un répo git pour le travail en groupe
- images/modif_parameters : graphes de $\log(v)=f(t)$ et $\log(v)=f(t)$ en ND pour des couples $(k,a)$ différent pour observer un comportement stable et instable (avec $\eta$ =1.0E-8), et graphe avec $\eta=1,0.10^{-11}$ (régime proche du quasi-statique).
- ajout d'une classe "dimensional parameters"
- portrait de phase $\dot{v} = f(v)$
- `rapport\theorie\theorie.tex` : obtention des 2 équations liant $\dot{v}$ et $\dot \theta$ à $v$ et $\theta$ que l'on utilise pour tracer $\log(v)=f(\theta)$ avec l'algorithme de Runge Kutta, détail des calculs écrit en LaTex.