#include <stdio.h>
#include <dlfcn.h>

void *cdlopen(char *name){
  void *handle;
  fprintf(stderr,"INFO: opening '%s'\n",name);
  handle = dlopen(name,RTLD_LAZY);
  if(! handle) {
    fprintf(stderr,"INFO: opening '%s' failed\n",name);
  }
  return(handle);
}