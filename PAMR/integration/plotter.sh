#! /bin/bash

for filename in mask*.dat; do

[ -e "$filename" ] || continue
number=`echo $filename | sed -e s/[^0-9]//g`

gnuplot<<END
set terminal png transparent enhanced font ',24' size 1200,2000 background rgb 'black'
set output "out.png"
set title "Grid $number" textcolor rgb 'white'

set xlabel "p" textcolor rgb 'white'
set ylabel "z" textcolor rgb 'white'
#set yrange [-20:20]
#set xrange [0:20]

set ytics 1
set xtics 1

set border lc rgb 'white'
set key off
set view map

unset colorbox

set palette model RGB defined (0 "white", 1 "black")
splot "mask$number.dat" u 1:2:3 with points pt 5 ps 1.5 palette

END

gnuplot<<END
set terminal png notransparent enhanced font ',24' size 1200,2000 background rgb 'black'
set output "out$number.png"
set title "Grid $number" textcolor rgb 'white'

set xlabel "p" textcolor rgb 'white'
set ylabel "z" textcolor rgb 'white'
#set yrange [-20:20]
#set xrange [0:20]

set ytics 1
set xtics 1

set border lc rgb 'white'
set key off
set view map

set palette model HSV
set palette rgb 3,2,2
splot "gfunc$number.dat" u 1:2:3 with points pt 5 ps 3 palette
END

convert out$number.png out.png -gravity center -composite out$number.png
rm out.png
xdg-open out$number.png

done
