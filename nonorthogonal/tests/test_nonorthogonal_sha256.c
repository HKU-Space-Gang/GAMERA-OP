#include "nonorthogonal_sha256.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  static const unsigned char expected[32] = {
      0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea,
      0x41, 0x41, 0x40, 0xde, 0x5d, 0xae, 0x22, 0x23,
      0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c,
      0xb4, 0x10, 0xff, 0x61, 0xf2, 0x00, 0x15, 0xad};
  gamera_no_sha256 context;
  unsigned char digest[32];
  gamera_no_sha256_init(&context);
  gamera_no_sha256_update(&context, "a", 1U);
  gamera_no_sha256_update(&context, "bc", 2U);
  gamera_no_sha256_final(&context, digest);
  if (memcmp(digest, expected, sizeof(expected)) != 0) {
    return 1;
  }
  puts("PASS SHA-256 incremental vector");
  return 0;
}
