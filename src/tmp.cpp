#include "std.h"
#include "mathfn.h"
#include "fileio.h"
#include "snp.h"



    struct Person
    {
        char name[50];
        int age;
        char phone[24];
    };

    int main()
    {
        int a = 9;
        char c = (char )a;
        c=c-48;
        cout << a << "," << c << endl;
    }

