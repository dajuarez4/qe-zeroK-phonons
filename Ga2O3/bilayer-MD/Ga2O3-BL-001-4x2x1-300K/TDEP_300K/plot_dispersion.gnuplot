set terminal pngcairo size 1100,700 enhanced font "Arial,13"
set output "Ga2O3-BL-001-4x2x1-TDEP-diagnostic.png"
set title "Ga2O3 (001) 4x2x1: diagnostic TDEP fit (2 MD frames)"
set ylabel "Frequency (THz)"
set xrange [0:0.261744]
set yrange [-35:35]
set xtics ("Γ" 0.0, "X" 0.084218, "S" 0.130872, "Y" 0.215090, "Γ" 0.261744)
set grid xtics ytics lc rgb "#d0d0d0"
set xzeroaxis linewidth 1.5 linecolor rgb "#202020"
unset key
plot for [column=2:31] "outfile.dispersion_relations" using 1:column \
    with lines lw 1.0 lc rgb "#285f9e"
