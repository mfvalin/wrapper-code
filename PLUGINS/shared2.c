#include <stdio.h>

char *entry_list[4] = { "Name_1","Name_2","Name_3",NULL} ;   // fully static method

int Name_1(int arg){
printf("Name_1: %d\n",arg);
return(arg);
}

int Name_2(int arg){
printf("Name_2: %d\n",arg);
return(arg);
}

int Name_3(int arg){
printf("Name_3: %d\n",arg);
return(arg);
}

