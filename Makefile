Csparce = -I/$(HOME)/software/CSparse/Include -L/$(HOME)/software/CSparse/Lib -lcsparse
GL = -L/usr/X11/lib -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -I/usr/X11/include



default :  timedep


timedep: timedep.c poisson1d_2.c steady.c
	$(CC) $(CLFAGS) -o $@ $^ $(Csparce) $(GL)

