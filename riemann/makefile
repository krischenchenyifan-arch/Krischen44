all: riemann main
	gcc main.o riemann.o -O3 -lm -o test.run
	gcc generate_data.o riemann.o -O3 -lm -o generate.run

main:
	gcc main.c -O3 -lm -c
	gcc generate_data.c -O3 -lm -c

riemann:
	gcc riemann.c -O3 -lm -c
