#!/usr/bin/env bash

echo "#### Beginning analysis pipeline ####"

echo "#### Cleaning directories: output, plots ####"
rm output/*
rm plots/*

echo "#### Running analysis scripts ####"
scripts/analyze_9nt_ss_ratios.py
scripts/analyze_11nt_ss_ratios.py
scripts/analyze_efficiency.py
scripts/analyze_junctions.py
scripts/analyze_libs.py
scripts/analyze_sums.py
scripts/compute_rnahybrid.py

echo "#### Running plotting scripts ####"
scripts/plot_correlations.py
scripts/plot_efficiency.py
scripts/plot_histograms.py
scripts/plot_junctions.py
scripts/plot_motifs.py
scripts/plot_scatter.py
scripts/plot_scatter_comput.py
scripts/plot_scatter_human.py
scripts/plot_validation.py
scripts/plot_venn.py
scripts/plot_barcodes.py
scripts/plot_clinical_variants.py
scripts/plot_regression.py

echo "#### Pipeline finished! ####"