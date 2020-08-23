#include <stdio.h>
#include <string.h>

char *entry_list[4] = { "name1","name2","name3",NULL} ;

int name1(int arg){
printf("C name1.1: %d\n",arg);
return(arg);
}

int name2(int arg){
printf("C name2.1: %d\n",arg);
return(arg);
}

int name3(int arg){
printf("C name3.1: %d\n",arg);
return(arg);
}

int get_symbol_number(){  // like fortran, function to get number of symbols, optional
  return(3);
}
