#include <udunits2.h>

main()
{
printf("\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\t%s=%d\n\n",
"integer, parameter :: UT_SUCCESS  ", UT_SUCCESS,          /* Success */
"integer, parameter :: UT_BAD_ARG  ", UT_BAD_ARG,         /* An argument violates the function's contract */
"integer, parameter :: UT_EXISTS  ", UT_EXISTS,          /* Unit, prefix, or identifier already exists */
"integer, parameter :: UT_NO_UNIT  ", UT_NO_UNIT,         /* No such unit exists */
"integer, parameter :: UT_OS  ", UT_OS,              /* Operating-system error.  See "errno". */
"integer, parameter :: UT_NOT_SAME_SYSTEM  ", UT_NOT_SAME_SYSTEM, /* The units belong to different unit-systems */
"integer, parameter :: UT_MEANINGLESS  ", UT_MEANINGLESS,     /* The operation on the unit(s) is meaningless */
"integer, parameter :: UT_NO_SECOND  ", UT_NO_SECOND,       /* The unit-system doesn't have a unit named "second" */
"integer, parameter :: UT_VISIT_ERROR  ", UT_VISIT_ERROR,     /* An error occurred while visiting a unit */
"integer, parameter :: UT_CANT_FORMAT  ", UT_CANT_FORMAT,     /* A unit can't be formatted in the desired manner */
"integer, parameter :: UT_SYNTAX  ", UT_SYNTAX,          /* string unit representation contains syntax error */
"integer, parameter :: UT_UNKNOWN  ", UT_UNKNOWN,         /* string unit representation contains unknown word */
"integer, parameter :: UT_OPEN_ARG  ", UT_OPEN_ARG,        /* Can't open argument-specified unit database */
"integer, parameter :: UT_OPEN_ENV  ", UT_OPEN_ENV,        /* Can't open environment-specified unit database */
"integer, parameter :: UT_OPEN_DEFAULT  ", UT_OPEN_DEFAULT,    /* Can't open installed, default, unit database */
"integer, parameter :: UT_PARSE  ", UT_PARSE  )          /* Error parsing unit specification */
;
}
