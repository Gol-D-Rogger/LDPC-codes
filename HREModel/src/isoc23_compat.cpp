#include <stdarg.h>
#include <stdio.h>

// Provide weak fallbacks for systems where glibc lacks the ISO C23 entry points.
extern "C" int __isoc23_sscanf(const char *str, const char *format, ...) __attribute__((weak));
extern "C" int __isoc23_fscanf(FILE *stream, const char *format, ...) __attribute__((weak));

extern "C" int __isoc23_sscanf(const char *str, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    int result = vsscanf(str, format, args);
    va_end(args);
    return result;
}

extern "C" int __isoc23_fscanf(FILE *stream, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    int result = vfscanf(stream, format, args);
    va_end(args);
    return result;
}