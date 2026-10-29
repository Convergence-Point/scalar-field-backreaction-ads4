# run from the repository root: gnuplot scripts/animation.gp
# columns of animation/frames.dat: x  Pi  Phi  A  delta
# frame times are in animation/times.dat (one per line, written by ./animation)
# loop -1: no loop extension in the GIF, so it plays once and stays on the last frame
# (loop 0 = repeat forever, loop N = repeat N more times)
#
# Optional settings, passed with gnuplot -e (all have defaults):
#   outfile='animation/wave.gif'   output file
#   maxframes=N                    render only the first N frames (default: all)
#   ymin=-1; ymax=1                fixed y range of the fields panel (default: from the data, see below)
#   gmin=-1.5; gmax=1.1            fixed y range of the metric panel (default: from the data)
#   weighted=1                     plot tan(x)*Pi and tan(x)*Phi instead of Pi and Phi. Away from the origin the
#                                  fields fall off like 1/tan(x), so with a fixed y range the unweighted pulse is
#                                  tiny except near x = 0; the weighted ones keep a roughly constant amplitude.
if (!exists("outfile"))   outfile = 'animation/wave.gif'
if (!exists("weighted"))  weighted = 0
wt(x) = (weighted ? tan(x) : 1.0)
set terminal gif animate delay 4 loop -1 size 900,700
set output outfile

# Values with |v| > cap (blow-up after a collapse) or NaN/Inf are ignored for the axis ranges
# and are not drawn, so one bad frame cannot make the y axis infinite. Raise `cap` if your
# fields legitimately exceed it.
cap = 1e3
ok(v) = (abs(v) < cap ? v : 1/0)

stats 'animation/frames.dat' using 1 nooutput
nframes = STATS_blocks - 1
if (exists("maxframes") && maxframes < nframes) { nframes = maxframes }

# Fields panel: ONE y range for the whole run, so the amplitude can be compared between frames.
# By default it is the range of Pi and Phi over the rendered frames, plus a margin of `pad`
# (fraction of the range) on both sides. Set ymin and ymax to override it.
pad = 0.08
if (!exists("ymin") || !exists("ymax")) {
    stats 'animation/frames.dat' index 0:(nframes-1) using (ok(wt($1)*$2)):(ok(wt($1)*$3)) nooutput
    flo = (STATS_min_x < STATS_min_y ? STATS_min_x : STATS_min_y)
    fhi = (STATS_max_x > STATS_max_y ? STATS_max_x : STATS_max_y)
    fm  = pad * (fhi - flo)
    if (fm < 1e-12) { fm = 0.1 * (abs(flo) + abs(fhi)) + 1e-12 }
    ymin = flo - fm
    ymax = fhi + fm
}
stats 'animation/frames.dat' index 0:(nframes-1) using (ok($4)) nooutput
amin = STATS_min
stats 'animation/frames.dat' index 0:(nframes-1) using (ok($5)) nooutput
dmin = STATS_min
dmax = STATS_max
# metric panel range (A <= 1; padded so that the flat A = 1, delta = 0 case is visible).
# Override with gmin and gmax; values outside the range are clipped at the panel edge.
if (!exists("gmin")) { gmin = (amin < dmin ? amin : dmin) - 0.05 }
if (!exists("gmax")) { gmax = (dmax > 1 ? dmax : 1) + 0.05 }

set xrange [0:pi/2]
# reference lines: vertical every pi/8 (labelled), horizontal major lines at the y tics plus 4 minor lines between
# them, and a darker line at zero.
set xtics ('0' 0, 'π/8' pi/8, 'π/4' pi/4, '3π/8' 3*pi/8, 'π/2' pi/2)
set mytics 5
set style line 100 linetype 0 linewidth 1 linecolor rgb '#7f7f7f'     # major grid lines (dotted, dark)
set style line 101 linetype 0 linewidth 1 linecolor rgb '#d0d0d0'     # minor grid lines (dotted, light)
set grid xtics ytics mytics linestyle 100, linestyle 101
set xzeroaxis linetype 1 linewidth 1 linecolor rgb '#555555'           # the line y = 0
set key top right

# Both panels keep one range for the whole run (fields: ymin..ymax above; metric: gmin..gmax above).

# frame times are read from animation/times.dat (one per line)
do for [k=0:nframes-1] {
    tk = word(system(sprintf('sed -n %dp animation/times.dat', k + 1)), 1) + 0
    set multiplot layout 2,1 title sprintf('AdS_4 radial wave,  t = %.3f', tk)

    set ylabel (weighted ? 'tan x · fields' : 'fields')
    unset xlabel
    set yrange [ymin:ymax]
    plot 'animation/frames.dat' index k using 1:(ok(wt($1)*$2)) with lines lw 2 title (weighted ? 'tan x · Π' : 'Π'), \
         'animation/frames.dat' index k using 1:(ok(wt($1)*$3)) with lines lw 2 title (weighted ? 'tan x · Φ' : 'Φ')

    set ylabel 'metric'
    set xlabel 'x  (r = L tan x)'
    set yrange [gmin:gmax]
    plot 'animation/frames.dat' index k using 1:(ok($4)) with lines lw 2 title 'A', \
         'animation/frames.dat' index k using 1:(ok($5)) with lines lw 2 title 'δ'

    unset multiplot
}
set output
