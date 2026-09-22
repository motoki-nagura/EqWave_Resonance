# GitHub publication manifest

This file lists the source files and directories that should be published.

## Rules

- Use absolute paths.
- One path per line.
- A file path selects that file.
- A directory path selects all files below that directory.
- Python files imported by selected Python files are added automatically when
  the imported module exists under `/A/data10/nagura/work/code`.
- Standard-library and third-party Python packages are not copied unless their
  source files are inside a configured local import root.
- Lines that do not begin with an absolute path are ignored by the script.
- Markdown bullet syntax and backticks are allowed.

## Publish

###  Figures
####      Figure 1
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/stats/spectrum/FIGURE_01`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/stats/spectrum/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/stats/spectrum/fig.pdf`

####      Figure 2
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/structures/FIGURE_02`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/structures/main.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/structures/func_draw_fig.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ObservedTaux/structures/fig.pdf`

####      Figure 3
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy/FIGURE_03`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy/fig.pdf`

####      Figures 4 and S2
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/FixedFrequency/results/box_bounded/FIGURE_04`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/FixedFrequency/results/box_bounded/FIGURE_S2`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/FixedFrequency/results/box_bounded/fig_Exp180dy+90dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/FixedFrequency/results/box_bounded/fig_Exp180dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/FixedFrequency/plot.py`

####      Figures 5, 12, and S3
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/FIGURE_05`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/FIGURE_12`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/FIGURE_S3`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/fig_Exp180dy+90dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/fig_ExpBottomDamp_180dy+90dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/results/box_bounded/fig_Exp180dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/k-m-omega/calc&plot.py`

####      Figure 6
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Filtering/filtering/Below500m/FIGURE_06`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Filtering/filtering/Below500m/fig_Exp180dy+90dy_T60dy.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Filtering/filtering/Below500m/main_filtering.py`

####      Figures 7 and 14
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/results/Exp_180dy+90dy/180-90-60day/x/FIGURE_07`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/results/Exp_180dy+90dy/180-90-60day/x/fig_3ptBox.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/results/Exp_180dy/180-180-90day/x/FIGURE_14`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/results/Exp_180dy/180-180-90day/x/fig_3ptBox.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/calc.py`

####      Figures 8, 9, and 10
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/FIGURE_08`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/FIGURE_09`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/FIGURE_10`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/fig_C1.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/fig_C2.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/Figure/fig_C3.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-90-60day/calc.py`

####      Figure 11
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy-Damp/FIGURE_11`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy-Damp/fig.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180+90dy-Damp/plot.py`

####      Figure 13
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180dy/FIGURE_13`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180dy/fig.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/omega/Averaged/Exp180dy/plot.py`

####      Figures 15, 16, and S4
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/FIGURE_15`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/FIGURE_16`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/FIGURE_S4`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/fig_C1.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/fig_C2.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/fig_C3.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/Figure/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_region/180-180-90day/calc.py`

####      Figure 17
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Mixed/FIGURE_17`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Mixed/calc.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Mixed/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Mixed/fig.pdf`

####      Figure A1
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Stratification/FIGURE_A1`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Stratification/fig.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Stratification/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Stratification/calc.sh`

####      Figure A2
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/KE+PE/FIGURE_A2`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/KE+PE/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/KineticEnergy/calc_parallel.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/KineticEnergy/merge.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/KineticEnergy/plot.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/PotentialEnergy/calc_parallel.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/PotentialEnergy/calc_InitialRho.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/PotentialEnergy/merge.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Energy/PotentialEnergy/plot.py`

####      Figure A3
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/wavelet/region/FIGURE_A3`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/wavelet/region/calc.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Spectrum/uvel/wavelet/region/plot.py`

####      Figure B1
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/significance/FIGURE_B1`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/significance/fig_3ptBox.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/plots/k3m3/k3m3_space/significance/plot.py`

####      Figure S1
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/plots/figure/tree/180-90-60dy/FIGURE_S1`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/plots/figure/tree/180-90-60dy/tree.pdf`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/plots/figure/tree/180-90-60dy/tree.tex`

####      Computation
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/config.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/main_solver.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/root_solver_utils.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-90-60dy/n1_1/params.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-90-60dy/n1_1/params_negative_m1.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-90-60dy/n1_1/params_near0_k1.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-180-90dy/n1_1/params.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-180-90dy/n1_1/params_negative_m1.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/Numerical_ver2/compute/params/N2_5e-6/180-180-90dy/n1_3/params.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/calc/main.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/calc_Synth/main.py`

####      Tools
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/CrossBispectrum.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/sub_clustering/clustering_2dArray.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/sub_find_index.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/sub_flip_omega_smth.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/sub_read_in_u_dudx_u.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/subs_k3m3_region_PlotTools.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/subs_k3m3_region_TargetRange.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/subs_lookup_table.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/subs_process_BS_ResTrio.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Bispectrum_tools/tools_k3m3_region.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/ResonantTrio_tools/sub_Read_ResonantTrios.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/Taux/sub_read_in_taux.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/calc_u_gradient.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/calc_w_interp.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/loop_uv.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/loop_w.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/merge.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/common/Bandpass_Frequencies.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/common/EqWaves_Rays.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/common/Fixed_Parameters.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/common/Overplot_DispersionRelation.py`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/common/Variance_from_PSD.py`
- `/S/data01/G4006/y0330/code/ofes/tools/post_post/ReadGrids/python/calc.py`
- `/A/data10/nagura/work/code/common/fft_spectrum_3d.py`
- `/A/data10/nagura/work/code/common/new_page.py`
- `/A/data10/nagura/work/code/common/regr_wgt.py`
- `/A/data10/nagura/work/code/common/spec_err.py`
- `/A/data10/nagura/work/code/common/sub_NetCDF_IO.py`
- `/A/data10/nagura/work/code/common/wavelet_subroutine.py`

###   Model

####      Initial Conditions
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Model/InitialConditions/1D/TS/ConstantN/calc.py`

####      Wind Forcing
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Model/Forcings/Steady/Idealized/add.sh`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Model/Forcings/Steady/Idealized/generate.jl`
- `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Model/Forcings/Steady/Idealized/loop.sh`

####      Grid generation
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/README`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/correct_kmt.f90`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/grids.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/ht.ctl`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/kmt.ctl`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/read_ht.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/read_kmt.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/read_kmt.out`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_correct_kmt`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_grids`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_read_ht`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_read_kmt`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_topog`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/run_topog_m`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/size.h`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/topog.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/topog_m.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/util.F`
- `/S/data01/G4006/y0330/code/ofes/topo/saved_settings/box/util_l100.F`

####      Grid data
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/grid.dta.out`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/kmt.ctl`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/kmt.dta`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/kmt.dta.bin`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/kmt.nc`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/results_grids`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/results_topog`
- `/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/size.h`

####     Model Parameters
- `/S/data01/G4006/y0330/code/ofes/scripts/box/box_bounded/run_ofes_ini.sh`
- `/S/data01/G4006/y0330/code/ofes/scripts/box/box_bounded/run_ofes_res.sh`

####     Compile Options
- `/S/data01/G4006/y0330/code/ofes/src/OFES/options/box/OPTIONS.OIFES_box_bounded`
- `/S/data01/G4006/y0330/code/ofes/src/OFES/options/box/OPTIONS.OIFES_box_bounded_BottomDrag`

<!-- Example: publish a file. -->
<!-- - `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools/some_file.py` -->

<!-- Example: publish an entire local directory. -->
<!-- - `/A/data10/nagura/work/code/themes/MidDepth/OFES_IdealExp/Tools` -->

<!-- Example: publish a file from /S/data01. -->
<!-- - `/S/data01/path/to/file.nc` -->
