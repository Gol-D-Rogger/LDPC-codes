#ifndef TEST_SVDPI_H
#define TEST_SVDPI_H

#include <stddef.h>

typedef struct test_sv_array_s {
  void *data;
  int len;
  size_t elem_size;
} test_sv_array_t;

typedef test_sv_array_t *svOpenArrayHandle;

static inline void *svGetArrayPtr(const svOpenArrayHandle handle) {
  return handle ? handle->data : NULL;
}

static inline void *svGetArrElemPtr(const svOpenArrayHandle handle, int index) {
  if (!handle || !handle->data || index < 0 || index >= handle->len)
    return NULL;
  return (void *)((char *)handle->data + ((size_t)index * handle->elem_size));
}

static inline int svHigh(const svOpenArrayHandle handle, int dimension) {
  (void)dimension;
  return (handle && handle->len > 0) ? (handle->len - 1) : -1;
}

#endif
