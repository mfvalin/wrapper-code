#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include <plugins.h>

int main() {
  fnptr fn;
  int arg=100;
  const charptr *table;
  void *p0, *p1, *p2, *p3, *p4, *p5;

  set_plugin_diag(1);

  p0 = load_plugin("libnoexist.so");
  if(p0) printf("pointer exists\n");

  p1 = load_plugin("libsharednone.so");
  if(p1) printf("pointer exists\n");

  p2 = load_plugin("libshared0.so");
  if(p2) printf("pointer exists\n");

  p3 = load_plugin("libshared1.so");
  if(p3) printf("pointer exists\n");

  p4 = load_plugin("libshared2.so");
  if(p4) printf("pointer exists\n");
  
  p5 = load_plugin("libsharedf1.so");  // cross test to fortran plugin
  if(p4) printf("pointer exists\n");
  
  printf(" function %s found at address %p\n","name1",fn=(fnptr)plugin_function(NULL,"name1"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","name2",fn=(fnptr)plugin_function(p3,"name2"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","name3",fn=(fnptr)plugin_function(p3,"name3"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","Name_1",fn=(fnptr)plugin_function(NULL,"Name_1"));
  if(fn) arg = (*fn)(arg) + 100;
  if(fn) printf(" function %s found at address %p\n","Name_1",fn=(fnptr)plugin_function(p4,"Name_2"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","Name_1",fn=(fnptr)plugin_function(p4,"Name_3"));
  if(fn) arg = (*fn)(arg) + 100;
  printf(" function %s found at address %p\n","name3f",fn=(fnptr)plugin_function(p5,"name3f"));
  if(fn) arg = (*fn)(&arg) + 100;
  printf(" function %s found at address %p\n","unadvertised",fn=(fnptr)plugin_function(p5,"unadvertised"));
  if(fn) arg = (*fn)(arg) + 100;

  table = plugin_function_names(p3);
  for(arg=0 ; arg<plugin_n_functions(p3) ; arg++) {
    printf(" fname libshared1.so = %s\n",table[arg]);
  }

  table = plugin_function_names(p4);
  for(arg=0 ; arg<plugin_n_functions(p4) ; arg++) {
    printf(" fname libshared2.so = %s\n",table[arg]);
  }

  table = plugin_function_names(p5);
  for(arg=0 ; arg<plugin_n_functions(p5) ; arg++) {
    printf(" fname libsharedf1.so = %s\n",table[arg]);
  }

  return(0);
}
