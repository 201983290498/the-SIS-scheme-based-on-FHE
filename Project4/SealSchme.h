#pragma once
#include<iostream>
#include "mymath.h"
#include"tools.h"
#include"SealSchme.h"
using namespace std;
using namespace seal;


/**
* print the parameters of BFV scheme;
*/
inline void print_parameters(std::shared_ptr<seal::SEALContext> context)
{
  // Verify parameters
  if (!context)
  {
    throw std::invalid_argument("context is not set");
  }
  auto& context_data = *context->key_context_data();

  /*
  Which scheme are we using?
  */
  std::string scheme_name;
  switch (context_data.parms().scheme())
  {
  case seal::scheme_type::BFV:
    scheme_name = "BFV";
    break;
  case seal::scheme_type::CKKS:
    scheme_name = "CKKS";
    break;
  default:
    throw std::invalid_argument("unsupported scheme");
  }
  std::cout << "/" << std::endl;
  std::cout << "| Encryption parameters :" << std::endl;
  std::cout << "|   scheme: " << scheme_name << std::endl;
  std::cout << "|   poly_modulus_degree: " <<
    context_data.parms().poly_modulus_degree() << std::endl;

  /*
  Print the size of the true (product) coefficient modulus.
  */
  std::cout << "|   coeff_modulus size: ";
  std::cout << context_data.total_coeff_modulus_bit_count() << " (";
  auto coeff_modulus = context_data.parms().coeff_modulus();
  std::size_t coeff_mod_count = coeff_modulus.size();
  for (std::size_t i = 0; i < coeff_mod_count - 1; i++)
  {
    std::cout << coeff_modulus[i].bit_count() << " + ";
  }
  std::cout << coeff_modulus.back().bit_count();
  std::cout << ") bits" << std::endl;

  /*
  For the BFV scheme print the plain_modulus parameter.
  */
  if (context_data.parms().scheme() == seal::scheme_type::BFV)
  {
    std::cout << "|   plain_modulus: " << context_data.
      parms().plain_modulus().value() << std::endl;
  }

  std::cout << "\\" << std::endl;
}


void getNoise(Decryptor& decryptor, Ciphertext xx) {
  cout << "Noise budget of x8x120_encrypted:" << decryptor.invariant_noise_budget(xx) << "bits" << endl;
}

Ciphertext getEnText(long long orgin, IntegerEncoder& encoder, Encryptor& encryptor) {
  Plaintext plain;
  Ciphertext en_text;
  encoder.encode(orgin, plain);
  encryptor.encrypt(plain, en_text);
  return en_text;
}

vector<Ciphertext> getEnText(vector<long long>orgin, IntegerEncoder& encoder,
  Encryptor& encryptor) {
  vector<Ciphertext> ans;
  Ciphertext en_text;
  for (int i = 0; i < orgin.size(); i++) {
    en_text = getEnText(orgin[i], encoder, encryptor);
    ans.push_back(en_text);
  }
  return ans;
}

void check(Decryptor& decryptor, Ciphertext xx) {
  int bits = decryptor.invariant_noise_budget(xx);
  if (bits > 0) {
    return;
  }
  else {
    cout << "the noise is out of range;\n";
    exit(0);
  }
}

int64_t getDeText(Decryptor& decryptor, IntegerEncoder& encoder, Ciphertext xx) {
  Plaintext plain;
  decryptor.decrypt(xx, plain);
  return encoder.decode_int64(plain);
}

void printDe(Decryptor& decryptor, IntegerEncoder& encoder, Ciphertext text) {
  Plaintext plain;
  decryptor.decrypt(text, plain);
  cout << encoder.decode_biguint(plain).to_dec_string() << endl;
}

// somewhate FHE scheme


vector<Ciphertext> getM(Evaluator& evaluator, RelinKeys& relinkeys, vector<Ciphertext> X,
  Ciphertext ONE, Decryptor& decryptor, Ciphertext& KT, IntegerEncoder& encoder) {//the decryptor is not necessary;
  vector<Ciphertext> M;
  for (int i = 0; i < X.size(); i++) {
    Ciphertext tem = ONE, tem1;
    for (int j = 0; j < X.size(); j++) {
      if (i == j)continue;
      evaluator.sub(X[i], X[j], tem1);
      evaluator.multiply_inplace(tem, tem1);
      evaluator.relinearize_inplace(tem, relinkeys);
      //check(decryptor, tem);
    }
    M.push_back(tem);
    evaluator.multiply_inplace(KT, tem);//KT may be vary large
    evaluator.relinearize_inplace(KT, relinkeys);
    check(decryptor, KT);
  }
  return M;
}

//the decryptor is not necessary;
vector<Ciphertext> getMCKKS(Evaluator& evaluator, RelinKeys& relinkeys, vector<Ciphertext> X,
  vector<Ciphertext> ONES, Decryptor& decryptor, CKKSEncoder& encoder, shared_ptr<seal::SEALContext>& context) {
  vector<Ciphertext> M;
  vector<Ciphertext> temVec;
  Ciphertext tem1;
  for (int i = 0; i < X.size(); i++) {
    temVec.clear();
    for (int j = 0; j < X.size(); j++) {
      if (i == j)continue;
      evaluator.sub(X[i], X[j], tem1);
      temVec.push_back(tem1);
    }
    tem1 = cal_mut(temVec, evaluator, relinkeys, context, ONES);
    M.push_back(tem1);
  }
  return M;
}

vector<Ciphertext> getMBFV(Evaluator& evaluator, RelinKeys& relinkeys, vector<Ciphertext> X,
  Ciphertext ONES, Decryptor& decryptor, BatchEncoder& encoder) {
  vector<Ciphertext> M;
  vector<Ciphertext> temVec;
  Ciphertext tem1,tem2;
  for (int i = 0; i < X.size(); i++) {
    temVec.clear();
    for (int j = 0; j < X.size(); j++) {
      if (i == j)continue;
      evaluator.sub(X[i], X[j], tem1);
      temVec.push_back(tem1);
    }
    tem2 = temVec[0];
    for (int i = 1; i < temVec.size(); i++){
      evaluator.multiply_inplace(tem2, temVec[i]);
      evaluator.relinearize_inplace(tem2, relinkeys);
    }
    //getNoise(decryptor, tem2);
    M.push_back(tem2);
  }
  return M;
}


vector<Ciphertext> getK(Evaluator& evaluator, RelinKeys& relinkeys, vector<Ciphertext> m,
  Ciphertext ONE, Decryptor& decryptor, Ciphertext P_se) {//the decryptor is not necessary;
  vector<Ciphertext> K;
  for (int i = 0; i < m.size(); i++) {
    Ciphertext tem = ONE;
    for (int j = 0; j < m.size(); j++) {
      if (i == j) continue;
      evaluator.multiply_inplace(tem, m[j]);
      evaluator.relinearize_inplace(tem, relinkeys);
    }
    //check(decryptor, tem);
    K.push_back(tem);
  }
  return K;
}

vector<Ciphertext> getACKKS(Evaluator& evaluator,
  RelinKeys& relinkeys,vector<Ciphertext>& Y,vector<Ciphertext>& m,
  vector<Ciphertext>& ONES, Decryptor& decryptor, shared_ptr<seal::SEALContext>& context) {//the decryptor is not necessary;
  vector<Ciphertext> A;
  vector<Ciphertext> temVec;
  Ciphertext tem;
  for (int i = 0; i < m.size(); i++) {
    temVec.clear();
    for (int j = 0; j < m.size(); j++) {
      if (i == j) continue;
      temVec.push_back(m[j]);
    }
    temVec.push_back(Y[i]);
    tem = cal_mut(temVec, evaluator, relinkeys, context, ONES);
    A.push_back(tem);
  }
  return A;
}

vector<Ciphertext> getABFV(Evaluator& evaluator,
  RelinKeys& relinkeys, vector<Ciphertext>& Y, vector<Ciphertext>& m,
  Ciphertext& ONES, Decryptor& decryptor) {//the decryptor is not necessary;
  vector<Ciphertext> A;
  vector<Ciphertext> temVec;
  Ciphertext tem;
  for (int i = 0; i < m.size(); i++) {
    temVec.clear();
    for (int j = 0; j < m.size(); j++) {
      if (i == j) continue;
      temVec.push_back(m[j]);
    }
    temVec.push_back(Y[i]);
    tem = temVec[0];
    for (int j = 1; j < temVec.size(); j++) {
      evaluator.multiply_inplace(tem, temVec[j]);
      evaluator.relinearize_inplace(tem, relinkeys);
    }
    A.push_back(tem);
  }
  return A;
}


/*void transform1(Decryptor& decryptor, IntegerEncoder& encoder, Ciphertext& k) {
  Plaintext plain;
  decryptor.decrypt(k, plain);
  BigUInt tem = encoder.decode_biguint(plain);
  cout << tem.to_dec_string() << endl;
  BigUInt ans;
  tem.divrem(251, ans);
  cout << ans.to_dec_string() << endl;
}*/

vector<int> recovery(vector<long long>& X, vector<long long>& Y, vector<long long>m,
  vector<long long>K, long long KT, long long invKT, size_t x, uint16_t plain_m, int k = 4) {

  cout << "--------------------------------------------\n---------------------------------------\n";
  cout << "Start recovery the pic by FHE;\n";
  cout << "Gennerate the scheme base information:\n";
  EncryptionParameters params(scheme_type::BFV);
  size_t poly_modulus_degree = x;
  params.set_poly_modulus_degree(poly_modulus_degree);
  cout << "the default xishu are: ";
  cout << CoeffModulus::BFVDefault(poly_modulus_degree).data() << endl;
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  params.set_plain_modulus(plain_m);
  auto context = SEALContext::Create(params);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  IntegerEncoder encoder(context);
  //encoder at local computer
  cout << "Start to encode the local information:\n";
  vector<Ciphertext> X_en = getEnText(X, encoder, encryptor);
  vector<Ciphertext> Y_en = getEnText(Y, encoder, encryptor);
  vector<Ciphertext> m_en = getEnText(m, encoder, encryptor);//the number of m is possible to be negative number
  vector<Ciphertext> K_en = getEnText(K, encoder, encryptor);
  Ciphertext KT_en = getEnText(KT, encoder, encryptor);
  Ciphertext invKT_en = getEnText(invKT, encoder, encryptor);


  cout << "upload X,Y,K,m,KT,invKT to the remote server:\n";

  vector<Ciphertext>A;//A = Ki*Yi
  Ciphertext tem;
  for (int i = 0; i < K_en.size(); i++) {
    evaluator.multiply(K_en[i], Y_en[i], tem);
    evaluator.relinearize_inplace(tem, relin_keys);
    check(decryptor, tem);
    A.push_back(tem);
  }
  Plaintext zeros = encoder.encode(0), one = encoder.encode(1);
  Ciphertext zeros_en, one_en;
  encryptor.encrypt(zeros, zeros_en);
  encryptor.encrypt(one, one_en);
  //mypow(evaluator, decryptor, encoder, KT_en, 250, one_en,relin_keys);
  vector<Ciphertext> D = getD(X_en, evaluator, decryptor, relin_keys, k, zeros_en, one_en);
  vector<Ciphertext> a;
  for (int i = 0; i < k; i++) {
    Ciphertext tem = zeros_en;
    Ciphertext tem1;
    for (int j = 0; j < k; j++) {
      evaluator.multiply(A[j], D[j * k + k - 1 - i], tem1);
      evaluator.add_inplace(tem, tem1);
    }
    evaluator.relinearize_inplace(tem, relin_keys);
    check(decryptor, tem);
    evaluator.multiply_inplace(tem, invKT_en);
    evaluator.relinearize_inplace(tem, relin_keys);
    check(decryptor, tem);
    a.push_back(tem);
  }

  cout << "download from  the center server;\n";
  cout << "\n\n\nstart to decoder locally\n";

  vector<int>ans_fl;
  long long tem_de;
  long long temP = 251;
  cout << "the orgin data of tem_de:\n";
  for (int i = 0; i < a.size(); i++) {
    tem_de = getDeText(decryptor, encoder, a[i]);
    cout << tem_de << " ";
    tem_de = ((tem_de % temP) % temP + temP) % temP;
    ans_fl.push_back(tem_de);
  }
  cout << "the recovery result of pixes using FHE:\n";
  for (int i = 0; i < ans_fl.size(); i++) {
    cout << ans_fl[i] << " ";
  }
  cout << endl;
  return ans_fl;
}



vector<Ciphertext> fullRecovery(vector<Ciphertext>& X_en, vector<Ciphertext>& Y_en,
  Ciphertext invKT_en, Encryptor& encryptor, Evaluator& evaluator, IntegerEncoder& encoder,
  Decryptor& decryptor, RelinKeys& relin_keys, Norm& norm, int k = 4) {

  /*
  vector<Ciphertext> m_en = getEnText(m, encoder, encryptor);//the number of m is possible to be negative number
  vector<Ciphertext> K_en = getEnText(K, encoder, encryptor);
  Ciphertext KT_en = getEnText(KT, encoder, encryptor);
  Ciphertext invKT_en = getEnText(invKT, encoder, encryptor);
  */

  Ciphertext zeros_en = norm.zeros_en, one_en = norm.one_en, KT_en, P_sen = norm.P_sen;
  KT_en = one_en;

  vector<Ciphertext>m_en = getM(evaluator, relin_keys, X_en, one_en, decryptor, KT_en, encoder);
  getNoise(decryptor, m_en[0]);
  vector<Ciphertext>K_en = getK(evaluator, relin_keys, m_en, one_en, decryptor, P_sen);
  getNoise(decryptor, K_en[0]);
  //cout << "the number KT is as follows:\n";
  //printDe(decryptor, encoder, KT_en);

  //decryptor.decrypt(KT_en, nowtem1);
  //cout << encoder.decode_int64(nowtem1) << endl;
 // encoder.encode(encoder.decode_int64(nowtem1), nowtem1);
  //encryptor.encrypt(nowtem1, nowtem);
  //invKT_en = mypow(evaluator,nowtem, decryptor, 249, one_en, relin_keys,encoder);
  //transform1(decryptor, encoder, invKT_en);


  vector<Ciphertext>A;//A = Ki*Yi
  Ciphertext tem;
  for (int i = 0; i < K_en.size(); i++) {
    evaluator.multiply(K_en[i], Y_en[i], tem);
    evaluator.relinearize_inplace(tem, relin_keys);
    check(decryptor, tem);
    // printDe(decryptor, encoder, tem);
    A.push_back(tem);
  }
  getNoise(decryptor, A[0]);
  vector<Ciphertext> D = getD(X_en, evaluator, decryptor, relin_keys, k, zeros_en, one_en);
  vector<Ciphertext> a;
  for (int i = 0; i < k; i++) {
    Ciphertext tem = zeros_en;
    Ciphertext tem1;
    for (int j = 0; j < k; j++) {
      evaluator.multiply(A[j], D[j * k + k - 1 - i], tem1);
      getNoise(decryptor, A[j]);
      getNoise(decryptor, D[j * k + k - 1 - i]);
      evaluator.add_inplace(tem, tem1);
      getNoise(decryptor, tem);
    }
    evaluator.relinearize_inplace(tem, relin_keys);
    check(decryptor, tem);

    // printDe(decryptor, encoder, tem);
     //here is  a problem that the type biguint is only judge the positive
    evaluator.multiply_inplace(tem, invKT_en);
    evaluator.relinearize_inplace(tem, relin_keys);
    //check(decryptor, tem);
    a.push_back(tem);
  }
  return a;
}


vector<Ciphertext> fullRecoveryCKKS(vector<Ciphertext>& X_en, vector<Ciphertext>& Y_en,
  Ciphertext invKT_en, Encryptor& encryptor, Evaluator& evaluator, CKKSEncoder& encoder,
  Decryptor& decryptor, RelinKeys& relin_keys, Norm& norm,Params& parms, shared_ptr<seal::SEALContext>& context) {
  int k = parms.getk();
  /*
  vector<Ciphertext> m_en = getEnText(m, encoder, encryptor);//the number of m is possible to be negative number
  vector<Ciphertext> K_en = getEnText(K, encoder, encryptor);
  Ciphertext KT_en = getEnText(KT, encoder, encryptor);
  Ciphertext invKT_en = getEnText(invKT, encoder, encryptor);
  */

  Ciphertext zeros_en = norm.ZERO,P_sen = norm.P_sen;

  vector<Ciphertext> ones = norm.ONES;
  vector<Ciphertext>m_en = getMCKKS(evaluator, relin_keys, X_en, ones, decryptor, encoder, context);
  //PrintDeVec(decryptor, encoder, m_en);
  vector<Ciphertext>A = getACKKS(evaluator, relin_keys, Y_en, m_en, ones, decryptor, context);
  //PrintDeVec(decryptor, encoder, A);

  vector<Ciphertext> D = getDCKKS(X_en, evaluator, decryptor, relin_keys, k, zeros_en,ones,context);
 /* PrintDeVec(decryptor, encoder, D);
  PrintDeVec(decryptor, encoder, invKT_en);
  print_line(3);*/
  vector<Ciphertext> a;
  for (int i = 0; i < k; i++) {
    Ciphertext tem = zeros_en;
    Ciphertext tem1;
    for (int j = 0; j < k; j++) {
      tem1 = cal_mut({ A[j], D[j * k + k - 1 - i] }, evaluator, relin_keys, context, ones);
      add(evaluator, tem, tem1,context);
    }
    //PrintDeVec(decryptor, encoder, tem);
    tem = cal_mut1({ tem, invKT_en }, evaluator, relin_keys, context, ones,decryptor,encoder);
    a.push_back(tem);
    //PrintDeVec(decryptor, encoder, tem);
  }
  return a;
}

vector<Ciphertext> fullRecoveryBFV(vector<Ciphertext>& X_en, vector<Ciphertext>& Y_en,
  Ciphertext invKT_en, Encryptor& encryptor, Evaluator& evaluator, BatchEncoder& encoder,
  Decryptor& decryptor, RelinKeys& relin_keys, Norm& norm, Params& parms, shared_ptr<seal::SEALContext>& context) {
  int k = parms.getk();
  /*
  vector<Ciphertext> m_en = getEnText(m, encoder, encryptor);//the number of m is possible to be negative number
  vector<Ciphertext> K_en = getEnText(K, encoder, encryptor);
  Ciphertext KT_en = getEnText(KT, encoder, encryptor);
  Ciphertext invKT_en = getEnText(invKT, encoder, encryptor);
  */

  Ciphertext zeros_en = norm.ZERO, P_sen = norm.P_sen;

  Ciphertext ones = norm.one_en;
  vector<Ciphertext>m_en = getMBFV(evaluator, relin_keys, X_en, ones, decryptor, encoder);
  //PrintDeVec(decryptor, encoder, m_en);
  getNoise(decryptor, m_en[0]);
  vector<Ciphertext>A = getABFV(evaluator, relin_keys, Y_en, m_en, ones, decryptor);
  //PrintDeVec(decryptor, encoder, A);
  getNoise(decryptor, A[0]);
  vector<Ciphertext> D = getDBFV(X_en, evaluator, decryptor, relin_keys, k, zeros_en, ones);
  /* PrintDeVec(decryptor, encoder, D);
   PrintDeVec(decryptor, encoder, invKT_en);
   print_line(3);*/
  vector<Ciphertext> a;
  for (int i = 0; i < k; i++) {
    Ciphertext tem = zeros_en;
    Ciphertext tem1;
    for (int j = 0; j < k; j++) {
      evaluator.multiply(A[j], D[j * k + k - 1 - i], tem1);
      evaluator.add_inplace(tem, tem1);
      getNoise(decryptor, tem);
    }
    evaluator.relinearize_inplace(tem, relin_keys);
    //PrintDeVec(decryptor, encoder, tem);
    evaluator.multiply_inplace(tem, invKT_en);
    evaluator.relinearize_inplace(tem, relin_keys);
    a.push_back(tem);
    //PrintDeVec(decryptor, encoder, tem);
  }
  return a;
}

// the improve scheme of CKKS where  the format of X is plantext;
vector<Ciphertext> fullRecoveryCKKS2(vector<double>& X, vector<Ciphertext>& Y_en,
  Ciphertext invKT_en, Encryptor& encryptor, Evaluator& evaluator, CKKSEncoder& encoder,
  Decryptor& decryptor, RelinKeys& relin_keys, Norm& norm, Params& parms, shared_ptr<seal::SEALContext>& context) {

  int k = parms.getk();
  /*
  vector<Ciphertext> m_en = getEnText(m, encoder, encryptor);//the number of m is possible to be negative number
  vector<Ciphertext> K_en = getEnText(K, encoder, encryptor);
  Ciphertext KT_en = getEnText(KT, encoder, encryptor);
  Ciphertext invKT_en = getEnText(invKT, encoder, encryptor);
  */

  Ciphertext zeros_en = norm.ZERO, P_sen = norm.P_sen;

  vector<Ciphertext> ones = norm.ONES;
  vector<double>m,D;
  {
    //solve the plaintext of M and D
    for (int i = 0; i < X.size(); i++) {
      double tem = 1;
      for (int j = 0; j < X.size(); j++) {
        if (i == j)continue;
        tem *= (X[i] - X[j]);
      }
      m.push_back(1.0/tem);
    }
    D = getD(X, k);
  }



  //PrintDeVec(decryptor, encoder, m_en);
  //vector<Ciphertext>A = getACKKS(evaluator, relin_keys, Y_en, m_en, ones, decryptor, context);
  //PrintDeVec(decryptor, encoder, A); the promotion of A is to delete division

  /* PrintDeVec(decryptor, encoder, D);
   PrintDeVec(decryptor, encoder, invKT_en);
   print_line(3);*/
  vector<Ciphertext> a;
  for (int i = 0; i < k; i++) {
    Ciphertext tem = zeros_en;
    Ciphertext tem1,tem_en;
    Plaintext tem_plain;
    for (int j = 0; j < k; j++) {
      double temA = m[j] * D[j * k + k - 1 - i];
      encoder.encode(temA, parms.getScale(), tem_plain);
      encryptor.encrypt(tem_plain, tem_en);
      tem1 = cal_mut({ Y_en[j],  tem_en}, evaluator, relin_keys, context, ones);
      add(evaluator, tem, tem1, context);
    }
    //PrintDeVec(decryptor, encoder, tem);
    //tem = cal_mut1({ tem, invKT_en }, evaluator, relin_keys, context, ones, decryptor, encoder);
    a.push_back(tem);
    //PrintDeVec(decryptor, encoder, tem);
  }
  return a;
}
