#include <stdio.h>
#include "qcbits.h"

static inline uint64_t RDTSC()
{
  unsigned int hi = 0, lo = 0;
  __asm__ volatile("rdtsc" : "=a" (lo), "=d" (hi));
  return ((uint64_t)hi << 32) | lo;
}


#ifdef OPENSSL
#include <openssl/evp.h>
void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  EVP_MD_CTX *hashctx = EVP_MD_CTX_new();
  uint8_t rnd[16];
  if(0 == _rdrand64_step ((unsigned long long *) rnd) ||
     0 == _rdrand64_step ((unsigned long long *) &(rnd[8]))){
    printf("Random Generation Failed\n");
    return;
  }

  if (!EVP_DigestInit_ex(hashctx, EVP_shake256(), NULL)
      || !EVP_DigestUpdate(hashctx, rnd, 16)
      || !EVP_DigestFinalXOF(hashctx, pointer, amount)){
        printf("Hash Calculation Failed\n");
        return;
      }
}

#else

void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  uint8_t rnd[8];
  int i = 0;
  while (i < amount) {
    if( i%8 == 0 && 0 == _rdrand64_step ((unsigned long long *) rnd)){
      printf("Random Generation Failed\n");
      continue;
    }
    pointer[i] = rnd[i%8];
    i++;
  }
}

#endif

#ifndef NUMBER_OF_TESTS
#define NUMBER_OF_TESTS 1000
#endif


int main(int argc, char const *argv[]) {
  word_t public_key[R_WORDS], symmetric_key[2][R_WORDS], ciphertext[R_WORDS], ciphertext2[R_WORDS],
         res_symmetric_key[R_WORDS];
  index_t private_key[2][W/2];
  uint64_t it_end[21] = {0};

  uint64_t start, end, keygen_time = 0, encrypt_time = 0, decrypt_time = 0, it;
  for (size_t ki = 0; ki < 25 + NUMBER_OF_TESTS; ki++) {

    /* Key Generation */
    start = RDTSC();
    if(1 == generate_keys(public_key, private_key)){
      printf("Key Generation Failure\n");
      continue;
    }
    end = RDTSC();
    if(ki >= 25) keygen_time += (end-start);

    /* Encryption */
    start = RDTSC();
    encrypt(public_key, symmetric_key, ciphertext);
    end = RDTSC();
    if(ki >= 25) encrypt_time += (end-start);

    /* Decryption */
    start = RDTSC();
    decrypt(private_key, ciphertext, res_symmetric_key);
    end = RDTSC();
    if(ki >= 25) decrypt_time += (end-start);

    /* Verification */
    for (size_t i = 0; i < 512; i++) {
      if(((uint64_t *)symmetric_key[0])[i] != ((uint64_t *)ciphertext)[i] ||
      ((uint64_t *)symmetric_key[1])[i] != ((uint64_t *)res_symmetric_key)[i]){
        printf("Decoding Failure %ld: %lx != %lx or %lx != %lx\n", i, ((uint64_t *)symmetric_key[0])[i], ((uint64_t *)ciphertext)[i],
        ((uint64_t *)symmetric_key[1])[i], ((uint64_t *)res_symmetric_key)[i]);
        break;
      }
    }
  }

  printf("Key Generation Time: %lf cycles\n", (double) keygen_time/NUMBER_OF_TESTS);
  printf("Encryption Time: %lf cycles\n", (double) encrypt_time/NUMBER_OF_TESTS);
  printf("Decryption Time: %lf cycles\n", (double) decrypt_time/NUMBER_OF_TESTS);

  return 0;
}
