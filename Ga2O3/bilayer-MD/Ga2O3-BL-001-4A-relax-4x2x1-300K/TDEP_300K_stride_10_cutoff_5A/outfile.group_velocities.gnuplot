 set terminal qt  size 500,350 enhanced font "CMU Serif,10"
 unset xtics
 set xtics ( "Γ" 0.0 ) 
set xtics add ("X"  0.084218  )
set xtics add ("S"  0.130872  )
set xtics add ("Y"  0.215090  )
set xtics add ("Γ"  0.261744  )
 set grid xtics lc rgb "#888888" lw 1 lt 0
 set xzeroaxis linewidth 0.1 linecolor 0 linetype 1
 set ytics scale 0.5
 set xtics scale 0.5
 set mytics 10
 unset key
 set ylabel "Group velocity (km/s)"
plot "outfile.group_velocities" u 1:2 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:3 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:4 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:5 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:6 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:7 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:8 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:9 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:10 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:11 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:12 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:13 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:14 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:15 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:16 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:17 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:18 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:19 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:20 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:21 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:22 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:23 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:24 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:25 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:26 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:27 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:28 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:29 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:30 w line lc rgb "#618712",\
 "outfile.group_velocities" u 1:31 w line lc rgb "#618712"
