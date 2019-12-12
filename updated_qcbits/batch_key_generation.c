

// TODO: Define N (the number of keys to generate)

// TODO: Insert this function in mainTest.c and call it from main
void testBatch(){
  word_t public_key[N][R_WORDS], symmetric_key[2][R_WORDS], ciphertext[R_WORDS], res_symmetric_key[R_WORDS];
  index_t private_key[N][2][W/2];
  uint64_t it_end[21] = {0};

  uint64_t start, end, keygen_time = 0, encrypt_time = 0, decrypt_time = 0, it;

  for (size_t ki = 0; ki < 25 + NUMBER_OF_TESTS; ki++) {
    if(!(ki %100)) {
      printf("%d\n", ki);
    }
    start = RDTSC();
    if(1 == batch_generate_keys(public_key, private_key)){
      printf("Key Generation Failure\n");
      continue;
    }
    end = RDTSC();
    if(ki >= 25) keygen_time += (end-start);
  }

  printf("Time: %lf\n", (double) keygen_time/NUMBER_OF_TESTS);
}

// TODO: Insert this function in qcbits.c and declare it in qcbits.h
int batch_generate_keys(word_t public_key[N][R_WORDS], index_t private_key[N][2][W/2]){

  word_t prod0[N][R_WORDS];
  word_t prod1[N][R_WORDS];
  word_t tmp[R_WORDS];
  // Generate the private key
  for(size_t i = 0; i < N; i++){
    generate_sparse_polynomial(W/2, private_key[i][0]);
    generate_sparse_polynomial(W/2, private_key[i][1]);
  }

  // PROD 0
  get_dense_polynomial_representation(private_key[0][0], W/2, prod0[0]);
  for(size_t i = 1; i < N; i++){
    multiply_sparse_polynomial(private_key[i][0], W/2, prod0[i - 1], prod0[i]);
  }

  // Inversion
  invert_polynomial(prod0[N - 1], prod1[N - 1]);

  // Inversion check
  multiply_polynomial(prod0[N - 1], prod1[N - 1], tmp);
  for (size_t i = 0; i < 511; i++) {
    if(((uint64_t *) tmp)[i]){
      printf("Erro %lx\n", ((uint64_t *) tmp)[i]);
      return 1;
    }
  }
  if(((uint64_t *) tmp)[511] != (((uint64_t) 1)<<63)){
    printf("*Erro %lx\n", ((uint64_t *) tmp)[511]);
    return 1;
  }

  // PROD 1
  for(size_t i = N - 2; i > 0; i--){
    multiply_sparse_polynomial(private_key[i + 1][0], W/2, prod1[i + 1], prod1[i]);
  }

  // Calculate Public keys
  for(size_t i = 0; i < N; i++){
    multiply_polynomial(prod1[i], prod0[i - 1], public_key[i]);
    multiply_sparse_polynomial(private_key[i][1], W/2, public_key[i], public_key[i]);
  }
  return 0;
}