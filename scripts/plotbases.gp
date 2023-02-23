# script: plotbases.gp
#
# plot stacked histogram from tbases.dat obtained with :
# 
# bioawk -c fastx '{print $comment}' \
#  V350134580_L01_e1_f91_demux/V350134580_L01_435301ABE103GEX_2.fq.gz \
#  | grep -v AAACCATTCTTCCCACCCTC \
#  |  /opt/biotools/SplitBarcode/dev/count_bases.pl
# SP@NC, 2023-02-22

# run with gnuplot plotbases.gp

set terminal jpeg giant font "Helvetica" 12

xdim=600
ydim=300
ar = ydim/xdim
set terminal png size xdim,ydim
set output 'test.png'
set size ratio ar

set output 'bases-substitutions.png'
set key outside
set key right top

set grid y
set style data histograms
set style histogram rowstacked
set boxwidth 0.5
set style fill solid 1.0 border -1
#set ytics 10 nomirror
#set yrange [0:130000]
set xlabel "Merged Barcode Position"
set ylabel "Base Frequencies"
#set ytics 10

plot 'tbases.dat' using 2 t "A", '' using 3 t "T", '' using 4 t "G", '' using 5 t "C", '' using 6:xtic(1) t "N"