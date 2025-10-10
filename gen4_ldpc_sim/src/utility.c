#include <time.h>
#include <stdio.h>
#include <string.h>

void print_time()
{
    time_t now = time(NULL);
    struct tm *tm = localtime(&now);
    printf("%s", asctime(tm));
}