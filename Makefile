all:
	mkdir -p out; cd validator; gcc -std=c11 -Wall -O2 ispd_validate.c -o ../out/ispd_validate -lm

