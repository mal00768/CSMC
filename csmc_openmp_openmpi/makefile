#make file for paringmp
# Mark Lingle
# 10/14/2013

OBJS = mainmp.o
CC = g++
DEBUG = -g
CFLAGS = -c $(DEBUG) -std=c++11 -O3 -fopenmp
LFLAGS = $(DEBUG) -I include -fopenmp

mainmp.x : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o mainmp.x

clean:
	\rm *.o *~ mainmp.x
