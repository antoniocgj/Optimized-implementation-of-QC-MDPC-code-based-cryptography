#include <x86intrin.h>
#include <stdint.h>

/* Security parameters */
#define R 32749
#define W 274
#define T 264

/* Implementation parameters */

#define word_t  __m512i
#define index_t uint16_t
#define R_WORDS 64

int generate_keys(word_t * public_key, index_t private_key[2][W/2]);
void encrypt(word_t * public_key, word_t dense_symmetric_key[2][R_WORDS], word_t * ciphertext);
int decrypt(index_t private_key[2][W/2], word_t * ciphertext, word_t * dense_symmetric_key);
extern void generate_random_bytes(uint64_t, uint8_t *);
