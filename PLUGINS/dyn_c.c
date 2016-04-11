#include <stdio.h>

int by_val(int arg) {
  fprintf(stderr,"by_val C arg=%d\n",arg);
  return arg;
}

int by_adr(int *arg) {
  fprintf(stderr,"by_adr C arg=%d\n",*arg);
  return *arg;
}
