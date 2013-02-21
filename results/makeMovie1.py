import sys,os,glob;
from scitools.std import movie;
if(len(sys.argv)<3):
    print "you forgot filename-base as commandline argument, and/or filetype";
else:
    filenamebase = sys.argv[1];
    filetype = sys.argv[2];
    movie(filenamebase + "*."+filetype, encoder='convert',fps=25, output_file='movie_'+filenamebase+'.gif');
    if(len(sys.argv)==4 and sys.argv[3]=='rm'):
        for i in glob.glob("%s*.%s"%(filenamebase,filetype)):
            os.remove(i);
