library(DiagrammeR)

grViz("
digraph SEM {

  graph [layout = dot, rankdir = LR]

  # -------------------------
  # EXOGENOUS VARIABLES
  # -------------------------
  Temp [label = 'Temperature Anomaly', shape = box]
  Precip [label = 'Antecedent Precipitation', shape = box]

  # -------------------------
  # NETWORK STATES (OBSERVED)
  # -------------------------
  Network0 [label = 'Baseline Network', shape = box]
  Network1 [label = 'Perturbed Network', shape = box]

  # -------------------------
  # DERIVED METRIC: RESISTANCE
  # -------------------------
  Resistance [label = 'Community Resistance\n(Δ Network Metrics)', shape = ellipse]

  # -------------------------
  # DYNAMICAL OUTCOME: ROBUSTNESS
  # -------------------------
  Robustness [label = 'Community Robustness\n(Topological and Stochastic)', shape = ellipse]

  # -------------------------
  # CAUSAL STRUCTURE
  # -------------------------

  Temp -> Network1
  Precip -> Network1

  Network0 -> Resistance
  Network1 -> Resistance

  Network1 -> Robustness

}
")