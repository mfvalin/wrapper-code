#include <stdio.h>
int test_main(int argc, char**argv)
{
  int i;
  for(i=0;i<=argc;i++) printf("arg %d = '%s'\n",i,argv[i]);
  return(0);
}

