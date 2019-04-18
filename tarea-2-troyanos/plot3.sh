#!\bin\bash

echo -e "set term \"png\"
set xrange [-$3:$3]
set yrange [-$3:$3]
set size square 
set output \"$1\"
plot \"$2\" u 1:2 t \"Sol\" w l, \"$2\" u 3:4 t \"Mercurio\" w l, \"$2\" u 5:6 w l" | gnuplot
