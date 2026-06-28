#ifndef MR_LRC_GF256_H
#define MR_LRC_GF256_H

#include <cstddef>
#include <stdint.h>
#include <vector>

namespace mrlrc {

uint8_t gf256_mul(uint8_t a, uint8_t b);
uint8_t gf256_pow(uint8_t a, unsigned int power);
uint8_t gf256_inv(uint8_t a);
bool gf256_invert_matrix(const std::vector<uint8_t> &in, std::vector<uint8_t> *out, int n);
bool gf256_matrix_is_invertible(const std::vector<uint8_t> &in, int n);
const char *gf256_backend();

}

#endif
