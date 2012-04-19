
GL = -L/usr/X11/lib -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -I/usr/X11/include

default : steady timedep

steady: steady.c
	$(CC) $(CLFAGS) -o $@ $^

timedep: timedep.c
	$(CC) $(CLFAGS) -o $@ $^ $(GL)

