#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <plugins.h>

typedef int (*fnptr)();             // pointer to function

int main() {
  fnptr fn;
  int arg=100;
  const char **table;
  void *p0, *p1, *p2, *p3, *p4, *p5;

  Set_plugin_diag(1);

  p0 = Load_plugin("libnoexist.so");
  if(p0) printf("pointer exists\n");

  p1 = Load_plugin("libsharednone.so");
  if(p1) printf("pointer exists\n");

  p2 = Load_plugin("libshared0.so");
  if(p2) printf("pointer exists\n");

  p3 = Load_plugin("libshared1.so");
  if(p3) printf("pointer exists\n");

  p4 = Load_plugin("libshared2.so");
  if(p4) printf("pointer exists\n");
  
  p5 = Load_plugin("libsharedf1.so");  // cross test to fortran plugin
  if(p4) printf("pointer exists\n");
  printf("\n");

  printf("blind lookup for name = name1\n");
  printf(" function %s found at address %p\n","name1",fn=(fnptr)Plugin_function(NULL,"name1"));
  if(fn) arg = (*fn)(arg) + 100;

  printf("targeted lookup in libshared1 for names = name2, name3\n");
  printf(" function %s found at address %p\n","name2",fn=(fnptr)Plugin_function(p3,"name2"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","name3",fn=(fnptr)Plugin_function(p3,"name3"));
  if(fn) arg = (*fn)(arg) + 100;

  printf("blind lookup for name = Name_1\n");
  printf(" function %s found at address %p\n","Name_1",fn=(fnptr)Plugin_function(NULL,"Name_1"));
  if(fn) arg = (*fn)(arg) + 100;

  printf("targeted lookup in libshared2 for names = Name_2, Name_3\n");
  if(fn) printf(" function %s found at address %p\n","Name_2",fn=(fnptr)Plugin_function(p4,"Name_2"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","Name_3",fn=(fnptr)Plugin_function(p4,"Name_3"));
  if(fn) arg = (*fn)(arg) + 100;

  printf("targeted lookup in libsharedf1 for names = name3f, unadvertised\n");
  printf(" function %s found at address %p\n","name3f",fn=(fnptr)Plugin_function(p5,"name3f"));
  if(fn) arg = (*fn)(&arg) + 100;
  printf(" function %s found at address %p\n","unadvertised",fn=(fnptr)Plugin_function(p5,"unadvertised"));
  if(fn) arg = (*fn)(arg) + 100;
  printf("\n");

  table = Plugin_function_names(p3);
  for(arg = 0 ; arg < Plugin_n_functions(p3) ; arg++) {
    printf(" plugin libshared1.so : %s\n",table[arg]);
  }
  printf("\n");

  table = Plugin_function_names(p4);
  for(arg = 0 ; arg < Plugin_n_functions(p4) ; arg++) {
    printf(" plugin libshared2.so : %s\n",table[arg]);
  }
  printf("\n");

  table = Plugin_function_names(p5);
  for(arg = 0 ; arg < Plugin_n_functions(p5) ; arg++) {
    printf(" plugin libsharedf1.so : %s\n",table[arg]);
  }
  printf("\n");

  return(0);
}
