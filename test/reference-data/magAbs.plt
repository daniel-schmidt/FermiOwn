reset
set xrange [0.1:0.3]
set yrange [0:1.1]
set arrow from 0.1864463,0 to 0.1864463,1.1 nohead
set xlabel "kappa"
set ylabel "absolute value of magnetization"
set key top left
plot 'magAbs8.dat' w lp t "8^3", 'magAbs10.dat' w lp t "10^3",'magAbs16.dat' w lp t "16^3", 'magAbs20.dat' w  lp t "20^3", '../magAbsMetropolis.dat' w lp t "Metropolis"
