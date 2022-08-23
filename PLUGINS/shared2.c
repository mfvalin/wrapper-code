#include <stdio.h>

char *EntryList_[4] = { "Name_1","Name_2","Name_3",NULL} ;   // fully static method

int Name_1(int arg){
printf("C Name_1.2: %d\n",arg);
return(arg);
}

int Name_2(int arg){
printf("C Name_2.2: %d\n",arg);
return(arg);
}

int Name_3(int arg){
printf("C Name_3.2: %d\n",arg);
return(arg);
}

void __attribute__ ((constructor)) Constructor2(void) {
   printf("plugin constructor for shared2\n");
}

void __attribute__ ((destructor)) Destructor2(void) {
   printf("plugin destructor for shared2\n");
}

