#include "qcbits.h"
/* Efficiency parameters */
#define DECODING_ITERATIONS 11
#define SUM_DEPTH 8


/* Polynomial Arithmetic Library */
#include "polyArithm.c"
// NOTE: This is not a generic library.

/* Internal use functions */
void generate_sparse_polynomial(int hammingWeight, index_t * indexes);
void decode(word_t[R_WORDS], index_t[2][W/2 + HWPA], word_t[R_WORDS],
            word_t[R_WORDS], word_t[SUM_DEPTH]);

/* Key Generation */
int generate_keys(word_t * public_key, index_t private_key_rel[2][W/2 + HWPA], index_t private_key_transposed[2][W/2]){

  index_t private_key[2][W/2];
  // Generate the private key
  generate_sparse_polynomial(W/2, private_key[0]);
  generate_sparse_polynomial(W/2, private_key[1]);

  // Transpose the private key. It will be used at the decryption.
  transpose_sparse_polynomial(private_key[0], W/2, private_key_transposed[0]);
  transpose_sparse_polynomial(private_key[1], W/2, private_key_transposed[1]);

  // Invert the private key
  word_t dense_private_key[2][R_WORDS];
  word_t inverted_private_key[R_WORDS];
  get_dense_polynomial_representation(private_key[0], W/2, dense_private_key[0]);
  get_dense_polynomial_representation(private_key[1], W/2, dense_private_key[1]);
  get_relative_sparse_representation(dense_private_key[0], private_key_rel[0], W/2 + HWPA);
  get_relative_sparse_representation(dense_private_key[1], private_key_rel[1], W/2 + HWPA);

  invert_polynomial(dense_private_key[0], inverted_private_key);

  // Check the inversion result
  multiply_sparse_polynomial(private_key_rel[0], W/2 + HWPA, inverted_private_key, public_key); 

  for (size_t i = 0; i < 511; i++) {
    if(((uint64_t *) public_key)[i]){
      printf("Erro %lx\n", ((uint64_t *) public_key)[i]);
      return 1;
    }
  }
  if(((uint64_t *) public_key)[511] != (((uint64_t) 1)<<63)){
    printf("*Erro %lx\n", ((uint64_t *) public_key)[511]);
    return 1;
  }

  // Calculate the public key
  multiply_sparse_polynomial(private_key_rel[1], W/2 + HWPA, inverted_private_key, public_key);
  return 0;
}


/* Symmetric key generation and encapsulation */
void encrypt(word_t * public_key, word_t dense_symmetric_key[2][R_WORDS], word_t * ciphertext){
  index_t symmetric_key[2][T/2];
  index_t symmetric_key_rel[T/2];

  // Generate two error vectors (symmetric key)
  generate_sparse_polynomial(T/2, symmetric_key[0]);
  generate_sparse_polynomial(T/2, symmetric_key[1]);

  get_dense_polynomial_representation(symmetric_key[0], T/2, dense_symmetric_key[0]);
  get_dense_polynomial_representation(symmetric_key[1], T/2, dense_symmetric_key[1]);

  get_relative_sparse_representation(dense_symmetric_key[1], symmetric_key_rel, T/2 + 31);

  // C = e[1] * G + e[0]
  multiply_sparse_polynomial(symmetric_key_rel, T/2 + 31, public_key, ciphertext);
  add_polynomial(ciphertext, dense_symmetric_key[0], ciphertext);
}


/* Symmetric key desencapsulation */
int decrypt(index_t private_key_rel[2][W/2 + HWPA], index_t private_key_transposed[2][W/2], word_t * ciphertext, word_t * dense_symmetric_key){
  index_t thr_var;
  uint64_t hw_s, final_it = DECODING_ITERATIONS;

  word_t partial_syndrome[R_WORDS];
  word_t syndrome[R_WORDS];
  word_t ciphertext_original[R_WORDS];

  index_t private_key_column[2][W/2 + HWPA];
  get_relative_representation(private_key_transposed[0], W/2, private_key_column[0], W/2 + HWPA);
  get_relative_representation(private_key_transposed[1], W/2, private_key_column[1], W/2 + HWPA);


  vector_set_value(dense_symmetric_key, 0);


  word_t thd[SUM_DEPTH];
  copy_polynomial(ciphertext, ciphertext_original);

  for (size_t i = 0; i < DECODING_ITERATIONS; i++) {
    // Syndrome Calculation
    multiply_sparse_polynomial(private_key_rel[0], W/2 + HWPA, ciphertext, partial_syndrome);
    multiply_sparse_polynomial(private_key_rel[1], W/2 + HWPA, dense_symmetric_key, syndrome);
    add_polynomial(partial_syndrome, syndrome, syndrome);

    hw_s = hamming_weight(syndrome);
    thr_var = ((int)(17.489 + 0.0043536 * (hw_s))) + 1;

#ifdef UNIFORM
    if(hw_s == 0){
      final_it = i;
      break;
    }
#endif

    // Decoding
    set_negative(thd, 1 + thr_var);
    decode(syndrome, private_key_column, ciphertext, dense_symmetric_key, thd);
  }
  add_polynomial(ciphertext, ciphertext_original, ciphertext); // isolating error
  return final_it;
}


void decode(word_t syndrome[R_WORDS], index_t private_key_column[2][W/2 + HWPA],
            word_t c[R_WORDS], word_t d[R_WORDS], word_t thr[SUM_DEPTH]){

  word_t rotated_syndrome[R_WORDS];
  word_t bit_count_vector[SUM_DEPTH][R_WORDS];
  __mmask8 is512;

  // Private Key first half
  for (size_t i = 0; i < SUM_DEPTH; i++) {
    vector_set_value(bit_count_vector[i], 0);
  }
  copy_polynomial(syndrome, rotated_syndrome);
  for (size_t i = 0; i < W/2 + HWPA; i++) {
    index_t idx = private_key_column[0][i];
    is512 = ((short)idx >> 15);
    multiply_monomial_polynomial((idx & ~is512) | (512 & is512), rotated_syndrome, rotated_syndrome);
    conditional_add_polynomial_to_bit_counter(rotated_syndrome, bit_count_vector, i, ~is512);
  }
  subtract_threshold_from_bit_counter(thr, bit_count_vector);
  addnot_polynomial(c, bit_count_vector[SUM_DEPTH - 1], c);

  // Private Key second half
  for (size_t i = 0; i < SUM_DEPTH; i++) {
    vector_set_value(bit_count_vector[i], 0);
  }
  copy_polynomial(syndrome, rotated_syndrome);
  for (size_t i = 0; i < W/2 + HWPA; i++) {
    index_t idx = private_key_column[1][i];
    is512 = ((short)idx >> 15);
    multiply_monomial_polynomial((idx & ~is512) | (512 & is512), rotated_syndrome, rotated_syndrome);
    conditional_add_polynomial_to_bit_counter(rotated_syndrome, bit_count_vector, i, ~is512);
  }
  subtract_threshold_from_bit_counter(thr, bit_count_vector);
  addnot_polynomial(d, bit_count_vector[SUM_DEPTH - 1], d);
}


void generate_sparse_polynomial(int hammingWeight, index_t * indexes){
  uint16_t repetition, gtR, i;
  do{
    // Generate indexes
    index_t indexes_e[W];
    generate_random_bytes(W * sizeof(index_t), (uint8_t *) indexes_e);

    // Check for repetitions
    i = 0;
    for (size_t k = 0; k < W; k++) {
      indexes_e[k] = indexes_e[k] & 0x7FFF;
      gtR = ((short)(R - indexes_e[k] -1) >> 15);
      repetition = 0;
      for (size_t j = 0; j < i; j++) {
        repetition |= ((indexes_e[k] ^ indexes[j]) - 1);
      }

      if((repetition|gtR) == 0xFFFF){ // repeated ou greater than R indexes
        continue;
      }

      indexes[i++] = indexes_e[k];
      if(i == hammingWeight){ // hamming Weight is public information
        return;
      }
    }
  }while (i < hammingWeight);
}