#include <stdio.h>

char *EntryList_[1] = {NULL} ;  // empty entry list

int name1(int arg){
printf("C name1.0: %d\n",arg);
return(arg);
}

int name2(int arg){
printf("C name2.0: %d\n",arg);
return(arg);
}

int name3(int arg){
printf("C name3.0: %d\n",arg);
return(arg);
}

void __attribute__ ((constructor)) Constructor0(void) {
   printf("plugin constructor for shared0\n");
}

void __attribute__ ((destructor)) Destructor0(void) {
   printf("plugin destructor for shared0\n");
}
