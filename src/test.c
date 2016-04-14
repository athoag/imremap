#include <stdio.h>
#include <string.h>

void dump_stringp(char *string, char *stringp)
{
}

int main()
{
    char string[] = "A B C";
    char *stringp = string;
    const char *delim = " \t";
    char *token;

    // stringp is updated to point to the next token 'B'
    token = strsep(&stringp, delim); 
    printf("token = '%s', ", token);
    if (stringp == NULL)
        printf("stringp == NULL\n");
    else
        printf("stringp - string = %d\n", stringp - string);

    // stringp is updated to point to the next token 'C'
    token = strsep(&stringp, delim); 
    printf("token = '%s', ", token);
    if (stringp == NULL)
        printf("stringp == NULL\n");
    else
        printf("stringp - string = %d\n", stringp - string);

    // In case no delimiter was found, stringp is made NULL.
    token = strsep(&stringp, delim); 
    printf("token = '%s', ", token);
    if (stringp == NULL)
        printf("stringp = NULL\n");
    else
        printf("stringp - string = %d\n", stringp - string);

    return 0;
}