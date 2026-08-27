#ifndef GAMERA_NONORTHOGONAL_SHA256_H
#define GAMERA_NONORTHOGONAL_SHA256_H

#include <stddef.h>
#include <stdint.h>

typedef struct gamera_no_sha256 {
  uint32_t state[8];
  uint64_t bytes;
  unsigned char block[64];
  size_t used;
} gamera_no_sha256;

void gamera_no_sha256_init(gamera_no_sha256 *context);
void gamera_no_sha256_update(gamera_no_sha256 *context, const void *data,
                             size_t bytes);
void gamera_no_sha256_final(gamera_no_sha256 *context,
                            unsigned char digest[32]);

#endif
