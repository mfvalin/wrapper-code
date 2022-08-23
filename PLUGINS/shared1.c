#include <stdio.h>
#include <string.h>

// list of names, NULL terminated, mandatory
char *EntryList_[4] = { "name1","name2","name3",NULL} ;

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

void __attribute__ ((constructor)) Constructor1(void) {
   printf("plugin constructor for shared1\n");
}

void __attribute__ ((destructor)) Destructor1(void) {
   printf("plugin destructor for shared1\n");
}
