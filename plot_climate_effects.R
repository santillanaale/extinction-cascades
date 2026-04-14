plot_climate_effects <- function(mods, climate_var, data) {
  
  climate_range <- seq(
    min(data[[climate_var]], na.rm = TRUE),
    max(data[[climate_var]], na.rm = TRUE),
    length.out = 100
  )
  
  preds <- lapply(names(mods), function(m) {
    predict_metric(mods[[m]], climate_var, climate_range) %>%
      mutate(metric = m)
  }) %>%
    bind_rows() %>%
    left_join(
      sig_table %>%
        filter(ClimateVar == climate_var),
      by = c("metric" = "Metric")
    )
  
  raw <- data %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c(
          # Pollinators
          "functional.complementarity.HL",
          "FunRedundancy.Pols",
          "mean.number.of.links.HL",
          "number.of.species.HL",
          
          # Plants
          "functional.complementarity.LL",
          "FunRedundancy.Plants",
          "mean.number.of.links.LL",
          "number.of.species.LL"
        )
      )
    )
  
  preds <- preds %>%
    mutate(metric = factor(metric, levels = levels(raw$metric)))
  
  ggplot(preds, aes(x = .data[[climate_var]], y = fit)) +
    geom_line(aes(linetype = signif), color = "black") +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "gray80", alpha = 0.4
    ) +
    geom_point(
      data = raw,
      aes(
        x = .data[[climate_var]],
        y = value,
        color = Site,
        shape = SampleRound
      ),
      inherit.aes = FALSE,
      alpha = 0.6
    ) +
    scale_linetype_manual(
      values = c(significant = "solid", ns = "dashed"),
      guide = "none"
    ) +
    facet_wrap(
      ~ metric,
      nrow = 2,
      scales = "free_y",
      labeller = labeller(metric = metric_labels)
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = climate_var,
      y = "Predicted community resistance metric",
      color = "Site"
    ) +
    guides(shape = "none")
}