/*  test program for object oriented c programming. This programm does:
    - create a testobject
    - assign a value to it
    - print that value
    - clone the object
    - compares original to new cloned object
    - creates new, different testobject
    - compares original to new, different object
*/
#include "cimplyobjects.h"
#include "testobject.h"
#include <petsc.h>

int main(){
void * obj1 = new(testobject, 5);
void * obj2 = clone(obj1);
void * obj3 = new(testobject,8);

printf("The value of object 1 is %i \nThe value of object 2 is %i \n The value of object 3 is %i \n ",getvalue(obj1), getvalue(obj2), getvalue(obj3));

if (differ(obj1,obj2)) printf("1 and 2 differ??\n");
if (differ(obj1, obj3)) printf("objects 1 and 3 differ.\n!");
delete(obj1);
delete(obj2);
delete(obj3);
return 0;
}
