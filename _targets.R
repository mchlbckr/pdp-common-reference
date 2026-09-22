library(targets)

tar_option_set(
  packages = c("stats", "ggplot2", "patchwork", "lpSolve"),
  format = "rds"
)

source(file.path("R", "support.R"))
source(file.path("R", "validity.R"))
source(file.path("R", "plot_theme.R"))
source(file.path("simulations", "s1_additive_correlation.R"))
source(file.path("simulations", "s4_off_support_perturbation.R"))
source(file.path("simulations", "s5_estimated_validity.R"))
source(file.path("simulations", "s2_nonlinear_manifold.R"))
source(file.path("simulations", "s3_interactions.R"))
source(file.path("simulations", "s6_baseline_comparison.R"))
source(file.path("simulations", "s7_mtcars_case_study.R"))
source(file.path("simulations", "s8_oracle_inference.R"))
source(file.path("simulations", "s9_estimated_validity_inference.R"))
source(file.path("simulations", "s10_outer_bootstrap.R"))
source(file.path("simulations", "s11_bootstrap_sample_size.R"))
source(file.path("simulations", "s12_validity_learner_robustness.R"))
source(file.path("simulations", "s13_airquality_case_study.R"))
source(file.path("simulations", "s14_experiment_matrix.R"))
source(file.path("simulations", "s15_boston_case_study.R"))
source(file.path("simulations", "s16_wine_quality_case_study.R"))
source(file.path("simulations", "s17_large_mixed_case_studies.R"))
source(file.path("simulations", "s18_relaxed_reference.R"))
source(file.path("simulations", "s19_calibration_convergence.R"))
source(file.path("simulations", "s20_alpha_sensitivity.R"))
source(file.path("simulations", "s21_dimension_sensitivity.R"))
source(file.path("simulations", "s22_rcpd_rate.R"))
source(file.path("simulations", "s23_composition_trap.R"))
source(file.path("simulations", "s24_constraint_comparison.R"))
source(file.path("simulations", "s25_error_decomposition.R"))
source(file.path("simulations", "s26_grid_resolution.R"))
source(file.path("simulations", "s27_application_uncertainty.R"))
source(file.path("simulations", "s28_full_split_stability.R"))
source(file.path("simulations", "s29_rcpd_comparison.R"))
source(file.path("simulations", "s30_lipschitz_wine.R"))
source(file.path("simulations", "s31_projection_generalization.R"))
source(file.path("simulations", "s32_bank_offsupport_stress.R"))
source(file.path("simulations", "s33_california_longitude_application.R"))
source(file.path("simulations", "s34_conditional_subgroup_comparison.R"))

list(
  tar_target(
    project_status,
    list(
      framework = "support diagnostics and common-reference feature effects",
      stage = "journal manuscript"
    )
  ),
  tar_target(
    simulation_s1,
    run_s1()
  ),
  tar_target(
    s1_feature_effects_figure,
    plot_s1(simulation_s1),
    format = "file"
  ),
  tar_target(
    s1_support_diagnostics_figure,
    plot_s1_support(simulation_s1),
    format = "file"
  ),
  tar_target(
    simulation_s4,
    run_s4()
  ),
  tar_target(
    simulation_s5,
    run_s5()
  ),
  tar_target(
    simulation_s2,
    run_s2()
  ),
  tar_target(
    simulation_s3,
    run_s3()
  ),
  tar_target(
    simulation_s6,
    run_s6()
  ),
  tar_target(
    case_study_s7,
    run_s7()
  ),
  tar_target(
    simulation_s8,
    run_s8()
  ),
  tar_target(
    simulation_s9,
    run_s9()
  ),
  tar_target(
    simulation_s10,
    run_s10()
  ),
  tar_target(
    simulation_s11,
    run_s11()
  ),
  tar_target(
    simulation_s12,
    run_s12()
  ),
  tar_target(
    case_study_s13,
    run_s13()
  ),
  tar_target(
    simulation_s14,
    run_s14()
  ),
  tar_target(
    s14_summary,
    summarize_s14(simulation_s14)
  ),
  tar_target(
    case_study_s15,
    run_s15()
  ),
  tar_target(
    case_study_s16,
    run_s16()
  ),
  tar_target(
    case_studies_s17,
    list(run_s17_bank(), run_s17_shoppers())
  ),
  tar_target(
    simulation_s18,
    run_s18()
  ),
  tar_target(
    s18_summary,
    summarize_s18(simulation_s18)
  ),
  tar_target(
    simulation_s19,
    run_s19()
  ),
  tar_target(
    s19_summary,
    summarize_s19(simulation_s19)
  ),
  tar_target(
    simulation_s20,
    run_s20()
  ),
  tar_target(
    s20_summary,
    summarize_s20(simulation_s20)
  ),
  tar_target(
    simulation_s21,
    run_s21()
  ),
  tar_target(
    s21_summary,
    summarize_s21(simulation_s21)
  ),
  tar_target(
    simulation_s22,
    run_s22()
  ),
  tar_target(
    s22_summary,
    summarize_s22(simulation_s22)
  ),
  tar_target(
    simulation_s23,
    run_s23()
  ),
  tar_target(
    s23_summary,
    summarize_s23(simulation_s23)
  ),
  tar_target(
    simulation_s24,
    run_s24()
  ),
  tar_target(
    s24_summary,
    summarize_s24(simulation_s24)
  ),
  tar_target(
    simulation_s25,
    run_s25()
  ),
  tar_target(
    s25_summary,
    summarize_s25(simulation_s25)
  ),
  tar_target(
    simulation_s26,
    run_s26()
  ),
  tar_target(
    s26_summary,
    summarize_s26(simulation_s26)
  ),
  tar_target(
    simulation_s27,
    run_s27(case_study_s16, case_studies_s17)
  ),
  tar_target(
    s27_summary,
    summarize_s27(simulation_s27)
  ),
  tar_target(
    simulation_s28,
    run_s28()
  ),
  tar_target(
    s28_summary,
    summarize_s28(simulation_s28)
  ),
  tar_target(
    simulation_s29,
    run_s29()
  ),
  tar_target(
    s29_summary,
    summarize_s29(simulation_s29)
  ),
  tar_target(
    simulation_s30,
    run_s30(case_study_s16)
  ),
  tar_target(
    s30_summary,
    summarize_s30(simulation_s30)
  ),
  tar_target(
    simulation_s31,
    run_s31()
  ),
  tar_target(
    s31_summary,
    summarize_s31(simulation_s31)
  ),
  tar_target(
    simulation_s32,
    run_s32(case_studies_s17[[1]])
  ),
  tar_target(
    s32_summary,
    summarize_s32(simulation_s32)
  ),
  tar_target(
    simulation_s33,
    run_s33()
  ),
  tar_target(
    s33_summary,
    summarize_s33(simulation_s33)
  ),
  tar_target(
    s33_gate,
    gate_s33(simulation_s33)
  ),
  tar_target(
    s33_audit,
    audit_s33(simulation_s33)
  ),
  tar_target(
    s33_reference_diagnostics,
    reference_diagnostics_s33(simulation_s33)
  ),
  tar_target(
    simulation_s34,
    run_s34()
  ),
  tar_target(
    s34_summary,
    summarize_s34(simulation_s34)
  ),
  tar_target(
    s3_interactions_figure,
    plot_s3(simulation_s3),
    format = "file"
  ),
  tar_target(
    s6_baseline_comparison_figure,
    plot_s6(simulation_s6),
    format = "file"
  ),
  tar_target(
    s7_mtcars_case_study_figure,
    plot_s7(case_study_s7),
    format = "file"
  ),
  tar_target(
    s8_oracle_inference_figure,
    plot_s8(simulation_s8),
    format = "file"
  ),
  tar_target(
    s9_estimated_validity_inference_figure,
    plot_s9(simulation_s9),
    format = "file"
  ),
  tar_target(
    s10_outer_bootstrap_figure,
    plot_s10(simulation_s10),
    format = "file"
  ),
  tar_target(
    s11_bootstrap_sample_size_figure,
    plot_s11(simulation_s11),
    format = "file"
  ),
  tar_target(
    s12_validity_learner_robustness_figure,
    plot_s12(simulation_s12),
    format = "file"
  ),
  tar_target(
    s13_airquality_case_study_figure,
    plot_s13(case_study_s13),
    format = "file"
  ),
  tar_target(
    s14_experiment_matrix_figure,
    plot_s14(s14_summary),
    format = "file"
  ),
  tar_target(
    s15_boston_case_study_figure,
    plot_s15(case_study_s15),
    format = "file"
  ),
  tar_target(
    s16_wine_quality_case_study_figure,
    plot_s16(case_study_s16),
    format = "file"
  ),
  tar_target(
    s17_large_mixed_case_studies_figure,
    plot_s17(case_studies_s17),
    format = "file"
  ),
  tar_target(
    s18_relaxed_reference_figure,
    plot_s18(simulation_s18),
    format = "file"
  ),
  tar_target(
    s19_calibration_convergence_figure,
    plot_s19(simulation_s19),
    format = "file"
  ),
  tar_target(
    s2_manifold_figure,
    plot_s2(simulation_s2),
    format = "file"
  ),
  tar_target(
    s4_off_support_figure,
    plot_s4(simulation_s4),
    format = "file"
  ),
  tar_target(
    s21_dimension_sensitivity_figure,
    plot_s21(s21_summary),
    format = "file"
  ),
  tar_target(
    s22_rcpd_rate_figure,
    plot_s22(s22_summary),
    format = "file"
  ),
  tar_target(
    s23_composition_trap_figure,
    plot_s23(simulation_s23),
    format = "file"
  ),
  tar_target(
    s24_constraint_comparison_figure,
    plot_s24(simulation_s24),
    format = "file"
  ),
  tar_target(
    s25_error_decomposition_figure,
    plot_s25(s25_summary),
    format = "file"
  ),
  tar_target(
    s26_grid_resolution_figure,
    plot_s26(s26_summary),
    format = "file"
  ),
  tar_target(
    s27_application_uncertainty_figure,
    plot_s27(simulation_s27),
    format = "file"
  ),
  tar_target(
    s28_full_split_stability_figure,
    plot_s28(simulation_s28),
    format = "file"
  ),
  tar_target(
    s29_rcpd_comparison_figure,
    plot_s29(s29_summary),
    format = "file"
  ),
  tar_target(
    s30_lipschitz_wine_figure,
    plot_s30(simulation_s30),
    format = "file"
  ),
  tar_target(
    s32_bank_offsupport_stress_figure,
    plot_s32(simulation_s32),
    format = "file"
  ),
  tar_target(
    s33_california_longitude_application_figure,
    plot_s33(simulation_s33),
    format = "file"
  ),
  tar_target(
    s33_california_reference_diagnostics_figure,
    plot_s33_reference_diagnostics(simulation_s33),
    format = "file"
  ),
  tar_target(
    s34_conditional_subgroup_comparison_figure,
    plot_s34(simulation_s34),
    format = "file"
  )
)
