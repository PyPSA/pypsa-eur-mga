# Geometric analysis of the near-optimal feasible space of energy system models

## Notes

https://pyomo.readthedocs.io/en/latest/working_models.html

quasi-monte-carlo integration with low-discrepancy sequences

https://stackoverflow.com/questions/42761370/set-initial-primal-and-dual-values-for-variables-pyomo

## Installation

Requires a `PyPSA-Eur` installation in a parallel folder.

```bash
../geometric-analysis % git clone https://github.com/pypsa/pypsa-eur.git
../geometric-analysis % git clone https://git.scc.kit.edu/FN/pypsa-eur-mga.git
```

## Description

Previous work on developing heuristics for transmission network expansion planning in low carbon energy system models motivated further analysis of the problem's solution space. The fact that different variants of proposed heuristics provided solutions very close to the optimal objective value, but differed significantly in their individual investment decisions (e.g. reinforcement of the transmission network or generator expansion), indicates that actually many near-optimal solutions exist. In fact, previous work has shown that many technologically diverse pathways for mitigating the climate crisis exist that keep costs within a pre-specified range [1]. 

Feasible, but sub-optimal solutions may be preferable for reasons that are difficult to quantify and, therefore, can not be contained in the model formulation [2]. Public acceptance of large infrastructure projects, risk, and ease of implementation within current regulatory frameworks are prime examples of considerations which are exogenous to most energy system models.

Instead of providing just a singular optimal solution, presenting multiple alternative solutions and pointing at features that persist across many near-optimal solutions can remedy the lack of certainty in energy system models (Hennen2017). Communicating model results in this way helps to identify must-haves (technology that is common to all near-optimal solutions) and must-avoids (technology that is not deployed in any near-optimal solution) [3].  

A common technique for determining multiple near-optimal solutions is called Modelling to Generate Alternatives (MGA) which uses the optimal solution as an anchor point to explore the surrounding decision space for maximally different solutions [2]. 

In this work I will apply a modified version of MGA to the PyPSA-Eur model of the European energy system. The objective function is encoded as a constraint, such that the feasible space is constrained by the optimal objective value plus some acceptable cost increase e. Subsequently, each investment variable (e.g. transmission lines, generators, and storage units) is both minimised and maximised observing the specified range (1 + e) of the objective function. This way, the whole near-optimal solution space is spanned. One can create further cuts for this hypercube by aggregating investment variables by type or location (e.g. minimize/maximize the sum of offshore wind generators in Germany).

One goal of this work is to find metrics that quantify the size of the near-optimal decision space. Additionally, it should be investigated how the extent of investment flexibility changes, when more ambitious emission reduction targets ranging from 80% to 100% are applied and various levels of slack ranging from 1% to 100% are permitted. Moreover, this work will also have to address methods that limit the computational effort of deriving many alternative solutions such as warmstart problem solving.

In summary, this work seeks to shed light on the geometry of the near-optimal feasible decision space of an energy system model with European scope.

### Literature

[1] Li, F. G. N., & Trutnevyte, E. (2017). Investment appraisal of cost-optimal and near-optimal pathways for the UK electricity sector transition to 2050. Applied Energy, 189, 89–109. https://doi.org/10.1016/j.apenergy.2016.12.047

[2] DeCarolis, J., Daly, H., Dodds, P., Keppo, I., Li, F., McDowall, W., … Zeyringer, M. (2017). Formalizing best practice for energy system optimization modelling. Applied Energy, 194, 184–198. https://doi.org/10.1016/j.apenergy.2017.03.001 

[3] Hennen, M., Lampe, M., Voll, P., & Bardow, A. (2017). SPREAD – Exploring the decision space in energy systems synthesis. Computers and Chemical Engineering, 106, 297–308. https://doi.org/10.1016/j.compchemeng.2017.06.002 

[4] DeCarolis, J. F. (2011). Using modeling to generate alternatives (MGA) to expand our thinking on energy futures. Energy Economics, 33(2), 145–152. https://doi.org/10.1016/j.eneco.2010.05.002 