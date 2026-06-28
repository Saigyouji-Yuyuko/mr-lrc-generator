#include "gf256.h"

#include <cstring>
#include <vector>

#ifdef MR_LRC_USE_ISAL
#include <erasure_code.h>
#endif

namespace mrlrc {

uint8_t gf256_mul(uint8_t a, uint8_t b)
{
#ifdef MR_LRC_USE_ISAL
    return gf_mul(a, b);
#else
    uint8_t product = 0;

    while (b != 0) {
        if ((b & 1U) != 0) {
            product ^= a;
        }

        uint8_t carry = static_cast<uint8_t>(a & 0x80U);
        a = static_cast<uint8_t>(a << 1U);
        if (carry != 0) {
            a ^= 0x1dU;
        }
        b = static_cast<uint8_t>(b >> 1U);
    }

    return product;
#endif
}

uint8_t gf256_pow(uint8_t a, unsigned int power)
{
    uint8_t result = 1;

    while (power != 0) {
        if ((power & 1U) != 0) {
            result = gf256_mul(result, a);
        }
        a = gf256_mul(a, a);
        power >>= 1U;
    }

    return result;
}

uint8_t gf256_inv(uint8_t a)
{
    if (a == 0) {
        return 0;
    }
    return gf256_pow(a, 254);
}

bool gf256_invert_matrix(const std::vector<uint8_t> &in, std::vector<uint8_t> *out, int n)
{
    if (n < 0 || out == nullptr) {
        return false;
    }
    if (static_cast<std::size_t>(n) * static_cast<std::size_t>(n) != in.size()) {
        return false;
    }
    out->assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0);
    if (n == 0) {
        return true;
    }

#ifdef MR_LRC_USE_ISAL
    std::vector<uint8_t> copy = in;
    return gf_invert_matrix(copy.data(), out->data(), n) == 0;
#else
    int width = n * 2;
    std::vector<uint8_t> work(static_cast<std::size_t>(n) * static_cast<std::size_t>(width), 0);

    for (int row = 0; row < n; row++) {
        std::memcpy(&work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width)],
                    &in[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)],
                    static_cast<std::size_t>(n));
        work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
             static_cast<std::size_t>(n) + static_cast<std::size_t>(row)] = 1;
    }

    for (int col = 0; col < n; col++) {
        int pivot = -1;
        for (int row = col; row < n; row++) {
            if (work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                     static_cast<std::size_t>(col)] != 0) {
                pivot = row;
                break;
            }
        }

        if (pivot < 0) {
            return false;
        }

        if (pivot != col) {
            for (int c = 0; c < width; c++) {
                std::swap(work[static_cast<std::size_t>(col) * static_cast<std::size_t>(width) +
                               static_cast<std::size_t>(c)],
                          work[static_cast<std::size_t>(pivot) * static_cast<std::size_t>(width) +
                               static_cast<std::size_t>(c)]);
            }
        }

        uint8_t inv_pivot =
            gf256_inv(work[static_cast<std::size_t>(col) * static_cast<std::size_t>(width) +
                           static_cast<std::size_t>(col)]);
        for (int c = 0; c < width; c++) {
            auto idx = static_cast<std::size_t>(col) * static_cast<std::size_t>(width) +
                       static_cast<std::size_t>(c);
            work[idx] = gf256_mul(work[idx], inv_pivot);
        }

        for (int row = 0; row < n; row++) {
            if (row == col) {
                continue;
            }

            uint8_t factor =
                work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                     static_cast<std::size_t>(col)];
            if (factor == 0) {
                continue;
            }

            for (int c = 0; c < width; c++) {
                work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                     static_cast<std::size_t>(c)] ^=
                    gf256_mul(factor,
                              work[static_cast<std::size_t>(col) * static_cast<std::size_t>(width) +
                                   static_cast<std::size_t>(c)]);
            }
        }
    }

    for (int row = 0; row < n; row++) {
        std::memcpy(&(*out)[static_cast<std::size_t>(row) * static_cast<std::size_t>(n)],
                    &work[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                          static_cast<std::size_t>(n)],
                    static_cast<std::size_t>(n));
    }

    return true;
#endif
}

bool gf256_matrix_is_invertible(const std::vector<uint8_t> &in, int n)
{
    std::vector<uint8_t> inverse;
    return gf256_invert_matrix(in, &inverse, n);
}

const char *gf256_backend()
{
#ifdef MR_LRC_USE_ISAL
    return "ISA-L";
#else
    return "scalar-fallback";
#endif
}

}
