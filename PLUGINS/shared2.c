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

