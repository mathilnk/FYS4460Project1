from scitools.std import *
import sys;
if(len(sys.argv)>2):
    base_name = sys.argv[1]
    num_files = eval(sys.argv[2])
    glued_name = base_name + "_glued.xyz"
    glued_file = open(glued_name, 'w');

    current = 0;

    for i in range(num_files):
        filename = base_name + str(current) + ".xyz";
        infile = open(filename, 'r');
        filestring = infile.read();
        glued_file.write(filestring);
        current+=1;
        infile.close();
    glued_file.close();
else:
    print "you forgot the commandline arguments.\nsys.argv[1] = basename\nsys.argv[2] = number of files"

