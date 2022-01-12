aprox: main.o splines.o points.o aproksymator_na_bazie.o libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o libge.a
	-rm *.o *.a

intrp: main.o splines.o points.o interpolator.o libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o libge.a
	-rm *.o *.a

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o
	-rm *.o

chebyshev: main.o splines.o points.o czebyszew.o libge.a
	$(CC) -o chebyshev main.o splines.o points.o czebyszew.o libge.a -lm
	-rm *.o *.a

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

czebyszew.o: 
	$(CC) -I gaus -c czebyszew.c

libge.a:
	cd gaus && make 

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
