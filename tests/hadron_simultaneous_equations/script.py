from subprocess import call
from os.path import isfile as does_exist

def create_script(i):
    index = str(i)
    with open('data/mass_' + index + '.dat') as f:
        mass = f.read()
    
    with open('gnuplot_s.gpi', 'w') as f:
        f.write('set term pngcairo size 1920, 480\n')
        f.write('unset key\n')
        f.write('set xrange[0:1.5]\n')
        f.write('set yrange[0:1.5]\n')
        f.write('set zrange[-1000:1000]\n')
        f.write('set output "frames/frame_' + index + '.png"\n')
        f.write('input_0 = "data/equation_0_' + index + '.dat"\n')
        f.write('input_1 = "data/equation_1_' + index + '.dat"\n')
        f.write('input_2 = "data/equation_2_' + index + '.dat"\n')
        f.write('curve_0 = "data/zero_reg_0_' + index + '.dat"\n')
        f.write('curve_1 = "data/zero_reg_1_' + index + '.dat"\n')
        f.write('curve_2 = "data/zero_reg_2_' + index + '.dat"\n')
        file_string = 'data/intersection_' + index + '.dat'
        if does_exist(file_string):
            f.write('intersec = "data/intersection_' + index + '.dat"\n')
#        f.write('set multiplot layout 1,3 title ' + mass + '\n')
#        f.write('splot input_0 w l, 0, curve_0 lc rgb "red", curve_1 lc rgb "yellow", curve_2 lc rgb "green"\n')
#        f.write('splot input_1 w l, 0, curve_0 lc rgb "red", curve_1 lc rgb "yellow", curve_2 lc rgb "green"\n')
#        f.write('splot input_2 w l, 0, curve_0 lc rgb "red", curve_1 lc rgb "yellow", curve_2 lc rgb "green"\n')
#        f.write('unset multiplot\n')
        f.write('set term pngcairo size 640, 480\n')
        f.write('set out "zeroes/zeroes_' + index + '.png\n')
        f.write('set title ' + mass + '\n')
        if does_exist(file_string):
            f.write('plot curve_0 lc rgb "red", curve_1 lc rgb "yellow", curve_2 lc rgb "green", intersec w points pt 7 lc rgb "black"\n')
        else:
            f.write('plot curve_0 lc rgb "red", curve_1 lc rgb "yellow", curve_2 lc rgb "green"')

            
def run():
    num_frames = 180
    
    i = 0
    while i < num_frames:
#        print('Frame: ' + str(i))
        create_script(i)
        call(["gnuplot", "gnuplot_s.gpi"])
        call(["rm", "gnuplot_s.gpi"])
        i = i + 1
        
    

if __name__ == '__main__':
    run()
