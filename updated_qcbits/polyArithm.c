#include <x86intrin.h>
#include <stdio.h>
#include <stdint.h>

#ifdef DEBUG
#define inline
#endif

/* Implementation types (Warning: POSIX types reserved style) */
#define word_t  __m512i
#define index_t uint16_t

/* Lengths of R*/
#ifdef R
  #if R != 32749
    #error Invalid R parameter. Polynomial Arithmetic Library built using (x^32749 - 1)
  #endif
#else
  #define R 32749
  #define R_WORDS 64
  #define SUM_DEPTH 8
#endif



/* Masks */

#define WORD_BITS (sizeof(word_t) * 8)
#define WORD_BITS_LOG 9
#define WORD_MOD 0x1FF
#define TAIL (R % WORD_BITS)
#define R_LOG 15

#define MOD_R(X) \
  (X & 0x7FFF) - (R & (((short)(R - (X & 0x7FFF) - 1)) >> 15))

/* Internal use functions */
static inline word_t bitShiftLeft512xmm (word_t, uint16_t);
static inline word_t bitShiftRight512xmm (word_t, uint16_t);
static inline word_t bitShiftLeft512xmmLT64 (word_t, word_t, word_t, word_t *);
static inline word_t bitShiftRight512xmmLT64 (word_t, word_t, word_t, word_t *);
static inline word_t bitShiftLeft512xmmC19 (word_t);
static inline word_t bitShiftRight512xmmC19 (word_t);
static inline word_t bitShiftLeft512xmmC493 (word_t);
static inline word_t bitShiftRight512xmmC493 (word_t);
static inline word_t bitShiftRight512xmmCarry (word_t data, index_t count, word_t * carryOut, word_t carryIn);
static inline word_t genMask (uint16_t);
static inline __mmask64 getBit(word_t, uint32_t);
static inline word_t conditionalSetBit(word_t, uint32_t, __mmask16);

#define compTailShiftLeft(X) bitShiftLeft512xmmC19(X)
#define compTailShiftRight(X) bitShiftRight512xmmC19(X)
#define tailShiftLeft(X) bitShiftLeft512xmmC493(X)
#define tailShiftRight(X) bitShiftRight512xmmC493(X)


void invert_polynomial(word_t A[R_WORDS], word_t B[R_WORDS]);
static inline void get_dense_polynomial_representation(index_t *sparse, int sparse_size, word_t *dense);
static inline void multiply_sparse_polynomial(index_t *sparse, int sparse_size, word_t *dense, word_t *out);
static inline void add_polynomial(word_t *in1, word_t *in2, word_t *out);
static inline void addnot_polynomial(word_t *in1, word_t *in2, word_t *out);
static inline void precompute_syndrome_rotations(word_t *in);
static inline void multiply_monomial_polynomial(index_t monomial, word_t *poly, word_t *out);
static inline void add_polynomial_to_bit_counter(word_t *poly, word_t bit_counter[SUM_DEPTH][R_WORDS], size_t it);
static inline void subtract_threshold_from_bit_counter(word_t thd[SUM_DEPTH], word_t bit_counter[SUM_DEPTH][R_WORDS]);
static inline void vector_set_value(word_t *vector, uint64_t value);
static inline void copy_polynomial(word_t * vector, word_t * vector_copy);
static inline void transpose_sparse_polynomial(index_t *sparse, int sparse_size, index_t *sparse_trans);
static inline uint64_t hamming_weight(word_t * poly);

void printVec(word_t * out, int size){
  for (size_t i = 0; i < size; i++) {
    printf("0x%lx, 0x%lx, 0x%lx, 0x%lx, 0x%lx, 0x%lx, 0x%lx, 0x%lx,",  ((uint64_t *)out)[8*i+0], ((uint64_t *)out)[8*i+1], ((uint64_t *)out)[8*i+2], ((uint64_t *)out)[8*i+3], ((uint64_t *)out)[8*i+4], ((uint64_t *)out)[8*i+5], ((uint64_t *)out)[8*i+6], ((uint64_t *)out)[8*i+7]);
  }
  printf("\n");
}

static inline uint64_t hamming_weight(word_t * poly){
  uint64_t total = 0;
  for (size_t i = 0; i < R_WORDS*8; i++) {
    total += _mm_popcnt_u64 (((uint64_t *)poly)[i]);
  }
  return total;
}

static inline void transpose_sparse_polynomial(index_t *sparse, int sparse_size, index_t *sparse_trans){
  for (size_t i = 0; i < sparse_size; i++) {
    index_t isZero = ((int)(sparse[i] - 1)) >> 31;
    sparse_trans[i] = (R - sparse[i]) & ~isZero;
  }
}


static inline void get_dense_polynomial_representation(index_t *sparse, int sparse_size,
                                         word_t *dense){

  vector_set_value(dense, 0);
  for (size_t i = 0; i < sparse_size; i++) {
    for (size_t j = 0; j < R_WORDS; j++) {
      __mmask16 indexMatch = ((int)((((sparse[i] >> WORD_BITS_LOG) ^ (R_WORDS - j - 1))) - 1)) >> 31;
      dense[j] = conditionalSetBit(dense[j], sparse[i] & WORD_MOD, indexMatch);
    }
  }
}

static inline void multiply_polynomial(word_t *p1, word_t *p2, word_t *out){
  #define R_WORDS_128 R_WORDS * 4
  word_t aux0_w[2*R_WORDS], aux1_w[2*R_WORDS], x;
  __m128i *p1_128, *p2_128;
  __m128i *aux0, *aux1;

  p1_128 = (__m128i *) p1;
  p2_128 = (__m128i *) p2;
  aux0 = (__m128i *) aux0_w;
  aux1 = (__m128i *) aux1_w;


  const __m128i zero = _mm_setzero_si128 ();
  for (size_t i = 0; i < 2*R_WORDS_128; i++) {
    aux0[i] = zero;
    aux1[i] = zero;
  }

  for (size_t i = 0; i < R_WORDS_128; i++) {
    for (size_t j = 0; j < R_WORDS_128; j++) {
      aux0[i+j  ] ^= _mm_clmulepi64_si128 (p1_128[i], p2_128[j], 0x00);
      aux1[i+j  ] ^= _mm_clmulepi64_si128 (p1_128[i], p2_128[j], 0x01);
      aux0[i+j+1] ^= _mm_clmulepi64_si128 (p1_128[i], p2_128[j], 0x11);
      aux1[i+j  ] ^= _mm_clmulepi64_si128 (p1_128[i], p2_128[j], 0x10);
    }
  }



  for (size_t i = 2*R_WORDS - 1; i > 0; i--) {
    x = _mm512_alignr_epi64 (aux1_w[i], aux1_w[i-1], 7);
    aux0_w[i] ^= x;
  }
  aux0_w[0] ^= _mm512_alignr_epi64 (aux1_w[0], _mm512_setzero_si512(), 7);


  // Reduction
  for (int i = R_WORDS - 1; i >= 0; i--) {
    x = tailShiftLeft(aux0_w[i + 1]) | compTailShiftRight(aux0_w[i]);
    out[i] = aux0_w[i + R_WORDS] ^ x;
  }
  out[0] &= genMask(WORD_BITS - TAIL);

  // INV(A) * INV(B) / x = INV(A*B)   --- INV = bit inversion
  word_t tmp, carryOut, one = _mm512_set1_epi64(1);
  carryOut = _mm512_alignr_epi64 (out[0], out[0], 7);
  for (size_t j = 1; j < R_WORDS; j++) {
    out[j] = bitShiftLeft512xmmLT64(out[j], one, carryOut, &carryOut);
  }
  out[0] = compTailShiftRight(out[0]);
  out[0] = bitShiftLeft512xmmLT64(out[0], one, carryOut, &carryOut);
  out[0] = compTailShiftLeft(out[0]);
  #undef R_WORDS_128
}

static inline void multiply_sparse_polynomial(index_t *sparse, int sparse_size, word_t *dense, word_t *out){
  word_t rotated_syndrome[R_WORDS];
  vector_set_value(out, 0);
  precompute_syndrome_rotations(dense);
  for (size_t i = 0; i < sparse_size; i++) {
    multiply_monomial_polynomial(sparse[i], dense, rotated_syndrome);
    add_polynomial(rotated_syndrome, out, out);
  }
}

static inline void add_polynomial(word_t *in1, word_t *in2, word_t *out){
  for (size_t i = 0; i < R_WORDS; i++) {
    out[i] = in1[i] ^ in2[i];
  }
}

static inline void addnot_polynomial(word_t *in1, word_t *in2, word_t *out){
  for (size_t i = 0; i < R_WORDS; i++) {
    out[i] = in1[i] ^ ~in2[i];
  }
}



#define ROTATE_WORDS(IN_VEC, OUT_VEC, MASK, AMOUNT) \
  OUT_VEC[0] = _mm512_mask_mov_epi64(IN_VEC[0], MASK, IN_VEC[AMOUNT] & tailMask); \
  for (j = 1; j < R_WORDS - AMOUNT; j++){\
    OUT_VEC[j] = _mm512_mask_mov_epi64(IN_VEC[j], MASK, IN_VEC[j + AMOUNT]);}\
  for (j = 0; j < AMOUNT; j++){\
    OUT_VEC[j + R_WORDS - AMOUNT] = _mm512_mask_or_epi64(IN_VEC[j + R_WORDS - AMOUNT], MASK, compTailShiftRight(IN_VEC[j]), tailShiftLeft(IN_VEC[j + 1]));}

word_t g_precomp[3][R_WORDS];

static inline void precompute_syndrome_rotations(word_t *in){
  const word_t tailMask = genMask(WORD_BITS - TAIL);
  index_t j;
  ROTATE_WORDS(in, g_precomp[0], 0xFF, 0x20);
}

static inline void multiply_monomial_polynomial(index_t monomial, word_t *poly, word_t *out){
  const word_t tailMask = genMask(WORD_BITS - TAIL);
  // Word rotation (permutating in the vector)
  index_t word_rotation_amount = monomial >> WORD_BITS_LOG;
  __mmask8 mask;
  index_t rpi, i, j;
  word_t tmp1[R_WORDS];

  i = R_LOG - 1;
  rpi = 1 << (i - WORD_BITS_LOG);
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  for (size_t ki = 0; ki < R_WORDS; ki++) {
    tmp1[ki] = _mm512_mask_mov_epi64(poly[ki], mask, g_precomp[0][ki]);
  }

  i = R_LOG - 2;
  rpi >>= 1;
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  ROTATE_WORDS(tmp1, out, mask, rpi);

  i = R_LOG - 3;
  rpi >>= 1;
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  ROTATE_WORDS(out, tmp1, mask, rpi);

  i = R_LOG - 4;
  rpi >>= 1;
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  ROTATE_WORDS(tmp1, out, mask, rpi);

  i = R_LOG - 5;
  rpi >>= 1;
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  ROTATE_WORDS(out, tmp1, mask, rpi);

  i = R_LOG - 6;
  rpi >>= 1;
  mask = -((word_rotation_amount & rpi) >> (i - WORD_BITS_LOG));
  ROTATE_WORDS(tmp1, out, mask, rpi);

  // Rotated inside the word
  index_t inside_rotation_amount = monomial & WORD_MOD;
  word_t outCarry, inCarry;
  word_t tmp = compTailShiftRight(out[0]) | tailShiftLeft(out[1]);
  bitShiftRight512xmmCarry(tmp, inside_rotation_amount, &outCarry, _mm512_setzero_si512());
  for (int z = R_WORDS - 1; z >= 0; z--) {
    inCarry = outCarry;
    out[z] = bitShiftRight512xmmCarry(out[z], inside_rotation_amount, &outCarry, inCarry);
  }
  out[0] &= tailMask;
}

static inline void add_polynomial_to_bit_counter(word_t *poly, word_t bit_counter[SUM_DEPTH][R_WORDS], size_t it){
  word_t carry[R_WORDS], carry_tmp;

  for (size_t i = 0; i < R_WORDS; i++) {
    carry[i] = poly[i] & bit_counter[0][i];
    bit_counter[0][i] ^= poly[i];
  }
  const int it_red[9] = {4, 3, 3, 2, 2, 1, 1, 0, 0};

  for (size_t j = 1; j < SUM_DEPTH - it_red[it>>4]; j++) {
    for (size_t i = 0; i < R_WORDS; i++) {
      carry_tmp = carry[i] & bit_counter[j][i];
      bit_counter[j][i] ^= carry[i];
      carry[i] = carry_tmp;
    }
  }

}

static inline void set_negative(word_t *v, uint16_t t)
{
  uint16_t t_neg = -t;
  const word_t C_N1 = _mm512_set1_epi64(-1);
  const word_t C_0 = _mm512_setzero_si512();

  for (int i = 0; i < SUM_DEPTH; i++){
    v[i] = _mm512_mask_mov_epi64(C_0,-((t_neg >> i) & 1), C_N1);
  }
}

static inline void subtract_threshold_from_bit_counter(word_t thd[SUM_DEPTH], word_t bit_counter[SUM_DEPTH][R_WORDS]){
  word_t carry[R_WORDS], carry_tmp;

  for (size_t i = 0; i < R_WORDS; i++) {
    carry[i] = thd[0] & bit_counter[0][i];
    bit_counter[0][i] = bit_counter[0][i] ^ thd[0];
  }

  for (size_t j = 1; j < SUM_DEPTH; j++) {
    for (size_t i = 0; i < R_WORDS; i++) {
      carry_tmp = (bit_counter[j][i] & thd[j]) | (carry[i] & (bit_counter[j][i] ^ thd[j]));
      bit_counter[j][i] = bit_counter[j][i] ^ thd[j] ^ carry[i];
      carry[i] = carry_tmp;
    }
  }
}

static inline void vector_set_value(word_t *vector, uint64_t value){
  for (size_t i = 0; i < R_WORDS; i++) {
    vector[i] = _mm512_set1_epi64(value);
  }
}

static inline void copy_polynomial(word_t * vector, word_t * vector_copy){
  for (size_t i = 0; i < R_WORDS; i++) {
    vector_copy[i] = vector[i];
  }
}

static inline __mmask64 getBit(word_t word, uint32_t bit){
  bit = WORD_BITS - 1 - bit; // the representation is inverted
  word_t cmp  = _mm512_set1_epi8(1 << (bit & 0x7));
  __mmask64 cmpRes = _mm512_mask_cmpeq_epi8_mask (1 << (bit >> 3), word & cmp, cmp);
  return -((cmpRes >> (bit >> 3)) & 1);
}

static inline word_t conditionalSetBit(word_t word, uint32_t bit, __mmask16 cond){
  bit = WORD_BITS - 1 - bit; // the representation is inverted
  word_t mask  = _mm512_set1_epi32(1 << (bit & 0x1F));
  return _mm512_mask_or_epi32 (word, (1 << (bit >> 5)) & cond, word, mask);
}


void invert_polynomial(word_t A[R_WORDS], word_t V[R_WORDS]) {
  uint32_t i, j, Fini;
  word_t nextOne;
  word_t delta = _mm512_set1_epi64(-1);

  const word_t C64 = _mm512_set1_epi64(64);
  const word_t C1  = _mm512_set1_epi64(1);
  const word_t CUBITS_TAIL = _mm512_set1_epi64(WORD_BITS - TAIL);
  const word_t zero = _mm512_setzero_si512();
  const word_t idx = _mm512_set1_epi64(7);
  const word_t one = _mm512_set_epi64(0, 0, 0, 0, 0, 0, 0, 1);

  word_t U[R_WORDS], U_old[R_WORDS];
  word_t F[R_WORDS], F_old[R_WORDS];
  word_t S[R_WORDS];
  word_t tmp, carryOut, innerCarry;
  word_t tailMask = genMask(WORD_BITS - TAIL);

  __mmask8 maskF0, maskDelta;

  // B = 1, C = 0, G = x^r - 1
  for (i = 0; i < R_WORDS; i++) {
    F[i] = A[i];
    U[i] = zero;
    V[i] = zero;
    S[i] = zero;
  }

  U[R_WORDS - 1] = bitShiftLeft512xmm(one, WORD_BITS - 1);
  S[R_WORDS - 1] = bitShiftLeft512xmm(one, WORD_BITS - 1);
  S[0] = bitShiftLeft512xmm(one, WORD_BITS - TAIL - 1); // x^r

  for (i = 0; i < (int) R*2 ; i++) {
    nextOne = _mm512_lzcnt_epi64 (F[R_WORDS - 1]);
    nextOne = _mm512_permutexvar_epi64 (idx, nextOne);
    Fini = 0; //i*0.97 / WORD_BITS;

    // F <- F/X // SHIFT RIGHT
    F[Fini] = bitShiftLeft512xmmLT64(F[Fini], nextOne, zero, &carryOut);
    for (j = Fini + 1; j < R_WORDS; j++) {
      F[j] = bitShiftLeft512xmmLT64(F[j], nextOne, carryOut, &carryOut);
    }

    carryOut = _mm512_alignr_epi64 (U[0], U[0], 7);
    for (j = 1; j < R_WORDS; j++) {
      U[j] = bitShiftLeft512xmmLT64(U[j], nextOne, carryOut, &carryOut);
    }
    U[0] = compTailShiftRight(U[0]);
    U[0] = bitShiftLeft512xmmLT64(U[0], nextOne, carryOut, &carryOut);
    U[0] = compTailShiftLeft(U[0]);

    delta -= nextOne;

    maskF0    = ~_mm512_cmpeq_epi64_mask (nextOne, C64);
    maskDelta = _mm512_mask_cmplt_epi64_mask (maskF0, delta, zero);

    for (j = Fini; j < R_WORDS; j++){
      F_old[j] = F[j];
      F[j] = _mm512_mask_xor_epi64 (F[j], maskF0, F[j], S[j]);
      S[j] = _mm512_mask_mov_epi64 (S[j], maskDelta, F_old[j]);
    }

    for (j = 0; j < R_WORDS; j++){
      U_old[j] = U[j];
      U[j] = _mm512_mask_xor_epi64 (U[j], maskF0, U[j], V[j]);
      V[j] = _mm512_mask_mov_epi64 (V[j], maskDelta, U_old[j]);
    }

    delta = _mm512_mask_sub_epi64 (delta, maskDelta, zero, delta);
  }
}

static inline word_t bitShiftLeft512xmmLT64 (word_t data, word_t count, word_t carryIn, word_t * carryOut){
  word_t innerCarry, out;
  const word_t C64 = _mm512_set1_epi64(64);
  // shift less than 64
  innerCarry = _mm512_alignr_epi64 (data, data, 7); // <<< 64
  *carryOut = innerCarry;
  innerCarry = _mm512_mask_mov_epi64 (innerCarry, 1, carryIn);
  innerCarry = _mm512_srlv_epi64 (innerCarry, C64 - count);
  out = _mm512_sllv_epi64 (data, count);
  out = _mm512_or_epi32 (out, innerCarry);
  return out;
}

static inline word_t bitShiftRight512xmmLT64 (word_t data, word_t count, word_t carryIn, word_t * carryOut){
  word_t innerCarry, out;
  const word_t C64 = _mm512_set1_epi64(64);
  // shift less than 64
  innerCarry = _mm512_alignr_epi64 (data, data, 1); // >>> 64
  *carryOut = innerCarry;
  innerCarry = _mm512_mask_mov_epi64 (innerCarry, 0x80, carryIn);
  innerCarry = _mm512_sllv_epi64 (innerCarry, C64 - count);
  out = _mm512_srlv_epi64 (data, count);
  out = _mm512_or_si512 (out, innerCarry);
  return out;
}

static inline word_t bitShiftRight512xmmCarry (word_t data, index_t count, word_t * carryOut, word_t carryIn){
  word_t innerCarry, out, countVet, idx, idx1;
  idx = _mm512_set_epi32(0xf, 0xe, 0xd, 0xc, 0xb, 0xa, 0x9, 0x8, 0x7, 0x6, 0x5, 0x4, 0x3, 0x2, 0x1, 0x0);

  countVet = _mm512_set1_epi8((count >> 5) & 0xE);
  idx1 = _mm512_add_epi32(idx, countVet);
  data = _mm512_permutexvar_epi32(idx1, data);
  *carryOut = data;
  data = _mm512_mask_blend_epi32(0xFFFF >> ((count >> 5) & 0xe), carryIn, data);
  // shift less than 64
  count = (count & 0x3F);
  innerCarry = _mm512_mask_blend_epi64(0xFE, carryIn, data);
  innerCarry = _mm512_alignr_epi64 (innerCarry, innerCarry, 1);
  innerCarry = _mm512_slli_epi64 (innerCarry, 64 - count);
  out = _mm512_srli_epi64 (data, count);
  out = _mm512_or_si512 (out, innerCarry);
  return out;
}


static inline word_t genMask (uint16_t n){
    word_t mask = _mm512_set1_epi8(-1);
    mask = bitShiftLeft512xmm(mask, n);
    return mask;
}

static inline word_t bitShiftLeft512xmm (word_t data, uint16_t count){
  word_t innerCarry, out, idx, countVet, idx1;
  idx = _mm512_set_epi32(0xf, 0xe, 0xd, 0xc, 0xb, 0xa, 0x9, 0x8, 0x7, 0x6, 0x5, 0x4, 0x3, 0x2, 0x1, 0x0);
  countVet = _mm512_set1_epi32((count >> 5) & 0xe);
  idx1 = _mm512_sub_epi32(idx, countVet);
  data = _mm512_maskz_permutexvar_epi32(0xFFFF << ((count >> 5) & 0xe), idx1, data);
  // shift less than 64
  count = (count & 0x3F);
  innerCarry = _mm512_maskz_permutexvar_epi32((uint16_t) (0xFFFF << 2), _mm512_sub_epi32(idx, _mm512_set1_epi32(2)), data); // << 64
  innerCarry = _mm512_srli_epi64 (innerCarry, 64 - count);
  out = _mm512_slli_epi64 (data, count);
  out = _mm512_or_si512 (out, innerCarry);
  return out;
}

static inline word_t bitShiftRight512xmm (word_t data, uint16_t count){
  word_t innerCarry, out, idx, countVet, idx1;
  idx = _mm512_set_epi32(0xf, 0xe, 0xd, 0xc, 0xb, 0xa, 0x9, 0x8, 0x7, 0x6, 0x5, 0x4, 0x3, 0x2, 0x1, 0x0);
  countVet = _mm512_set1_epi32((count >> 5) & 0xe);
  idx1 = _mm512_add_epi32(idx, countVet);
  data = _mm512_maskz_permutexvar_epi32(0xFFFF >> ((count >> 5) & 0xe), idx1, data);
  // shift less than 64
  count = (count & 0x3F);
  innerCarry = _mm512_maskz_permutexvar_epi32(0xFFFF >> 2, _mm512_add_epi32(idx, _mm512_set1_epi32(2)), data); // >> 64
  innerCarry = _mm512_slli_epi64 (innerCarry, 64 - count);
  out = _mm512_srli_epi64 (data, count);
  out = _mm512_or_si512 (out, innerCarry);
  return out;
}


static inline word_t bitShiftLeft512xmmC19 (word_t data){
  word_t innerCarry, out;
  const word_t C00 = _mm512_setzero_si512();
  innerCarry = _mm512_alignr_epi64 (data, C00, 7); // <<< 64
  innerCarry = _mm512_srli_epi64 (innerCarry, 45);
  out = _mm512_slli_epi64 (data, 19);
  out = _mm512_or_epi32 (out, innerCarry);
  return out;
}

static inline word_t bitShiftRight512xmmC19 (word_t data){
  word_t innerCarry, out;
  const word_t C00 = _mm512_setzero_si512();
  innerCarry = _mm512_alignr_epi64 (C00, data, 1); // >>> 64
  innerCarry = _mm512_slli_epi64 (innerCarry, 45);
  out = _mm512_srli_epi64 (data, 19);
  out = _mm512_or_si512 (out, innerCarry);
  return out;
}

static inline word_t bitShiftLeft512xmmC493 (word_t data){
  const word_t C00 = _mm512_setzero_si512();
  data = _mm512_alignr_epi64 (data, C00, 1);
  return _mm512_slli_epi64 (data, 45);
}

static inline word_t bitShiftRight512xmmC493 (word_t data){
  const word_t C00 = _mm512_setzero_si512();
  data = _mm512_alignr_epi64 (C00, data, 7);
  return _mm512_srli_epi64 (data, 45);
}
