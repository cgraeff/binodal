#!/bin/bash

cd data/
list=`echo */`

for dir in $list; do

    cd $dir
    
    echo "   [$dir]"
    
    echo "set term pngcairo size 1024,768" > gnuplot_script.gpi

    echo "set yrange[-200:1400]" >> gnuplot_script.gpi

    echo "set title \"$dir\"" >> gnuplot_script.gpi
    echo "set xlabel \"barionic chemical potential\"" >> gnuplot_script.gpi
    echo "set ylabel \"pressure\"" >> gnuplot_script.gpi

    echo "do for [i=0:59] {" >> gnuplot_script.gpi

    echo "hadron_pressure_file=sprintf(\"hadron_pressure_%d.dat\",i)" >> gnuplot_script.gpi
    echo "quark_pressure_file=sprintf(\"quark_pressure_%d.dat\",i)" >> gnuplot_script.gpi
    echo "binodal_point_file=sprintf(\"binodal_point_file_%d.png\",i)" >> gnuplot_script.gpi
    echo "set out binodal_point_file" >> gnuplot_script.gpi
    echo "plot hadron_pressure_file w l t \"hadron\", quark_pressure_file w l t \"quark\", \"\" u 1:3 w l t \"mu_I\"" >> gnuplot_script.gpi
    
    echo "point_mass_file=sprintf(\"point_mass_file_%d.png\",i)" >> gnuplot_script.gpi
    echo "set out point_mass_file" >> gnuplot_script.gpi
    echo "plot hadron_pressure_file u 1:4 w lp pt 6 t \"hadron mass\"" >> gnuplot_script.gpi
    
    echo "}" >> gnuplot_script.gpi
        
    echo "set yrange[-1:2]" >> gnuplot_script.gpi
    
    echo "do for [i=0:59] {" >> gnuplot_script.gpi
    
    echo "hadron_pressure_file=sprintf(\"hadron_pressure_%d.dat\",i)" >> gnuplot_script.gpi
    
    echo "point_dens_file=sprintf(\"point_dens_file_%d.png\",i)" >> gnuplot_script.gpi
    echo "set out point_dens_file" >> gnuplot_script.gpi
    echo "plot hadron_pressure_file u 1:5 w lp pt 6 t \"proton dens\", \"\" u 1:6 w lp pt 6 t \"neutron dens\"" >> gnuplot_script.gpi

    echo "}" >> gnuplot_script.gpi
    
    gnuplot gnuplot_script.gpi
    
    rm gnuplot_script.gpi
    
    if [ ! -e "../../graph/$dir/" ];
    then
        mkdir -p "../../graph/$dir/";
    fi
    
    mv *.png "../../graph/$dir/"
    
    cd ..
done
