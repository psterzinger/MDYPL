plot_with_insets <- function(estimates, pt, type, cols = c("#CF4446", "#00AD9A", "#FB9A06"),
                             alpha_fac = 500, p_size = 0.5, vp_h = 0.14, vp_w = vp_h * 1.1) {
    kappa_gamma <- unique(estimates[c("kappa", "gamma")])
    insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        ep <- estimates |>
            subset(kappa == ckappa & gamma == cgamma)
        p_alpha <- alpha_fac * (1 - ckappa) / nrow(ep)
        ep <- ep |>
            plot_results(p_alpha = p_alpha, p_size = p_size, type = type, cols = cols) +
            geom_hline(aes(yintercept = 0), linetype = 3, col = "grey")
        insets_estimates[[wh]] <- tibble(x = ckappa,
                                         y = cgamma,
                                         plot = list(ep))
    }
    out <- plot_pt(pt, max_kappa = 1.0) +
        geom_point(data = kappa_gamma, aes(kappa, gamma),
                   pch = 18, size = 1,
                   col = "grey")
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        out <- out +  geom_plot(data = insets_estimates[[wh]],
                                aes(x = x, y = y, label = plot),
                                vp.width = vp_w,
                                vp.height = vp_h,
                                vjust = "bottom", hjust = "left")
    }
    out
}

## For a single kappa/gamma combination
plot_results <- function(summaries, cols = c("#CF4446", "#00AD9A", "#FB9A06"),
                         p_alpha = 0.2, p_size = 1,
                         type = c("estimate_vs_truth", "estimate_and_truth")) {
    type <- match.arg(type)
    lims <- with(summaries, range(c(truth, estimate)))
    if (type == "estimate_vs_truth") {
        p1 <- ggplot(summaries) +
            geom_point(aes(x = truth, y = estimate), alpha = p_alpha, size = p_size, col = cols[3]) +
            geom_abline(aes(intercept = 0, slope = 1), col = "black", lwd = 0.5) +
            geom_smooth(aes(x = truth, y = estimate),
                        method = "lm", formula = "y ~ x",
                        se = FALSE, col = cols[2], lwd = 0.5) +
            coord_cartesian(x = lims, y = lims)
    }
    if (type == "estimate_and_truth") {
        sd <- summaries |> group_by(kappa, gamma, truth) |>
            summarize(e_mean = mean(estimate)) |>
            inner_join(summaries, c("kappa", "gamma", "truth"))
        p1 <- ggplot(summaries) +
            geom_point(aes(x = parameter, y = estimate), alpha = p_alpha, size = p_size, col = cols[3]) +
            geom_line(aes(x = parameter, y = truth), col = "black", lwd = 0.5) +
            geom_line(data = sd, aes(x = parameter, y = e_mean), col = cols[2], lwd = 0.5) +
            coord_cartesian(y = lims)
    }
    p1 + theme_minimal() +
        theme(legend.position = "none") +
        theme(
            plot.background = element_rect(fill= 'transparent', color = "grey"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_blank()
        ) +
        labs(x = NULL, y = NULL)
}



plot_pt <- function(PT, max_kappa = NULL, cols = c("#F1F1F1", "#FFFFFF")) {
    PT <- PT[order(PT$kappa), ]
    max_gamma <- max(PT$gamma)
    max_kappa <- ifelse(is.null(max_kappa), max(PT$kappa), max_kappa)
    polygon1 <- data.frame(kappa = c(PT$kappa, max_kappa, max_kappa),
                           gamma = c(PT$gamma, 0, max_gamma))
    polygon2 <- data.frame(kappa = c(0, 0, PT$kappa),
                           gamma = c(0, max_gamma, PT$gamma))
    ggplot(PT) +
        geom_polygon(data = polygon1, aes(kappa, gamma), fill = cols[1]) +
        geom_polygon(data = polygon2, aes(kappa, gamma), fill = cols[2]) +
        geom_line(aes(kappa, gamma), col = "grey") +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = expression(kappa), y = expression(gamma))
}
