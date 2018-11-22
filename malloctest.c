#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

int main(){
	int *mem;
	unsigned long mem_size;
	mem = malloc(sizeof(unsigned int) * 6103515625);
	mem_size = malloc_usable_size(mem);
	printf("sizeof(mem) = %lld\n", mem_size);
}
