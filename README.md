# MIG-seismi-modelisation

_Objectif : modéliser la sismicité induite par injection de fluide en géothermie en Alsace._

## Jour 1
- création d'un répo git pour le travail en groupe
- images/modif_parameters : graphes de $\log(v)=f(t)$ et $\log(v)=f(t)$ en ND pour des couples $(k,a)$ différent pour observer un comportement stable et instable (avec $\eta$ =1,0.10^{-8}), et graphe avec $\eta=1,0.10^{-11}$ (régime proche du quasi-statique).
- ajout d'une classe "dimensional parameters"
- portrait de phase $\dot{v} = f(v)$, avec échelle de couleur pour la progression du temps
- `rapport\theorie\theorie.tex` : obtention des 2 équations liant $\dot{v}$ et $\dot \theta$ à $v$ et $\theta$ que l'on utilise pour tracer $\log(v)=f(\theta)$ avec l'algorithme de Runge Kutta, détail des calculs écrit en LaTex.
- _Segall, P., & Lu, S. (2015). Injection-induced seismicity : Poroelastic and earthquake nucleation effects. Journal of Geophysical Research: Solid Earth, 120(7), 5082‑5103. https://doi.org/10.1002/2015JB012060_
 : (cf. page 3) permet de modéliser le puit d'injection avec un nouveau paramètre (la distance r de la zone d'intérêt de la faille avec la zone d'injection). Intéressant car les effets thermo/poro élastiques sont des causes de micro-sismicité en géothermie profonde (cf p134 guide des bonnes pratiques).
