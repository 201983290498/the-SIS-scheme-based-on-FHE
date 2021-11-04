#pragma once
#include "baseClass.h"
#include"SealSchme.h"

using namespace std;
using namespace seal;

void getShareByHE(size_t x, uint16_t plain_m, Picture& source,
  vector<SharePic>& shares, Params& parms, Result& finalRes, bool printAns = true) {
  cout << "Start to generate secret shares by FHE(BFV):" << endl;
  vector<Ciphertext> X_en;
  EncryptionParameters params(scheme_type::BFV);
  size_t poly_modulus_degree = x;
  params.set_poly_modulus_degree(poly_modulus_degree);
  cout << "the default xishu are: ";
  cout << CoeffModulus::BFVDefault(poly_modulus_degree).data() << endl;
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  params.set_plain_modulus(plain_m);
  auto context = SEALContext::Create(params);
  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  IntegerEncoder encoder(context);

  vector<ll> y = source.getSec();
  //encrype the originPic
  vector<Plaintext> pix_plain;

  cout << "the process is encoding orginal Pictures:" << endl;
  double dur;
  clock_t st, end;
  st = clock();
  Plaintext tem1;
  Ciphertext tem2;
  for (int i = 0; i < y.size(); i++) {
    encoder.encode(y[i], tem1);
    encryptor.encrypt(tem1, tem2);
    pix_plain.push_back(tem1);
    source.pushCipher(tem2);
  }

  vector<Ciphertext> pix_en = source.getSecEn();

  //encoder and encryted the index of X
  for (int i = 0; i < shares.size(); i++) {
    vector<ll> xid = shares[i].getX();
    for (int j = 0; j < xid.size(); j++) {
      encoder.encode(xid[j], tem1);
      encryptor.encrypt(tem1, tem2);
      shares[i].x_en.push_back(tem2);
    }
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.encode0 = (dur / CLOCKS_PER_SEC);
  cout << "the encoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;

  cout << "the Process-generate Shares by FHE is beginning:" << endl;
  st = clock();
  vector <Ciphertext>temy;
  for (int index = 0; index * parms.getk() < y.size(); index = index + 1) {
    temy.clear();
    //get Every group of YCilpherText; 
    for (int i = 0; i < parms.getk(); i++) {
      temy.push_back(pix_en[index * parms.getk() + i]);
    }
    for (int pn = 0; pn < shares.size(); pn++) {
      Ciphertext& Tem = shares[pn].addNewPixByCipher(temy, evaluator, relin_keys,index);
      //getNoise(decryptor, Tem);
    }
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.generate0 = (dur / CLOCKS_PER_SEC);
  cout << "the Getting Shares' Process has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;


  cout << "the state of decoding the shares is beginning:" << endl;
  st = clock();
  //show the image
  for (int i = 0; i < shares.size(); i++) {
    if(printAns)
      cout << "the " << i << "th shares' pixes are:";
    shares[i].showShare(decryptor, encoder, printAns);
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.decode0 = (dur / CLOCKS_PER_SEC);
  cout << "the state of decoding the shares has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
}

Picture recoryShare(Params& picParms, vector<SharePic> uploadShares,
  DecTools& tools, Result& finalRes, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  cout << "--------------------------------------------\n---------------------------------------\n";
  cout << "Start Fully recovery the pic by FHE;\n";
  cout << "Gennerate the scheme base information:\n";
  EncryptionParameters params(scheme_type::BFV);
  size_t poly_modulus_degree = degree;
  params.set_poly_modulus_degree(poly_modulus_degree);
  cout << "the default xishu are: ";
  cout << CoeffModulus::BFVDefault(poly_modulus_degree).data() << endl;
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  params.set_plain_modulus(p_mo);
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
  cout << "Start to encode the local shares:\n";
  cout << "waiting...." << endl;
  double dur;
  clock_t start, end;
  start = clock();
  for (int i = 0; i < uploadShares.size(); i++) {
    for (int j = 0; j < uploadShares[i].X.size(); j++) {
      encryptor.encrypt(encoder.encode(uploadShares[i].X[j]), uploadShares[i].x_en[j]);
      encryptor.encrypt(encoder.encode(uploadShares[i].fx[j]), uploadShares[i].fx_en[j]);
    }
  }
  end = clock();
  dur = (double)(end - start);
  cout << "encoding procedure is finished.";
  cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
  finalRes.encode1 = (dur / CLOCKS_PER_SEC);

  Norm norm(encoder, encryptor, picParms);
  vector<Ciphertext>X, Ys;
  Ciphertext invKT;
  Picture recon_pic(picParms);
  cout << "Start reconstructing Pics procedure:" << endl;
  start = clock();
  for (int group = 0; group < uploadShares[0].X.size(); group++) {
    encryptor.encrypt(encoder.encode(tools.invKT[group]), invKT);
    X.clear(), Ys.clear();
    for (int i = 0; i < uploadShares.size(); i++) {
      X.push_back(uploadShares[i].x_en[group]);
      Ys.push_back(uploadShares[i].fx_en[group]);
    }
    vector<Ciphertext> a = fullRecovery(X, Ys, invKT, encryptor, evaluator, encoder, decryptor, relin_keys, norm, picParms.getk());
    for (int i = 0; i < a.size(); i++) {
      recon_pic.pushCipher(a[i]);
    }
  }
  end = clock();
  dur = (double)(end - start);
  cout << "reconstructing Pics procedure is finished." << endl << endl;
  cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl;
  finalRes.generate1 = (dur / CLOCKS_PER_SEC);

  recon_pic.setPic(recon_pic.DecryPic(decryptor, encoder, evaluator, norm, finalRes, printAns));
  return recon_pic;
}

/// <summary>
/// CKKS schemes of getShare
/// </summary>
/// <param name="x"></param>
/// <param name="plain_m"></param>
/// <param name="source"></param>
/// <param name="shares"></param>
/// <param name="parms"></param>
/// <param name="finalRes"></param>
void getShareByCKKS(size_t x, Picture& source, vector<SharePic>& shares,
  Params& parms, Result& finalRes,bool printAns = true) {
  cout << "Start to generate secret shares by FHE(BFV):" << endl;
  vector<Ciphertext> X_en;
  EncryptionParameters params(scheme_type::CKKS);
  size_t poly_modulus_degree = x;
  params.set_poly_modulus_degree(poly_modulus_degree);
  vector<int>coeff_modulus;

  coeff_modulus.push_back(min(parms.getBaseLenOfModulus() + 20,60));
  for (int i = 0; i < parms.getMaxLevel(); i++) {
    coeff_modulus.push_back(parms.getBaseLenOfModulus());
  }
  coeff_modulus.push_back(min(parms.getBaseLenOfModulus() + 20,60));

  params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree,coeff_modulus));

  auto context = SEALContext::Create(params);
  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  CKKSEncoder encoder(context);

  vector<double> y = source.getSecFromDouble();
  vector<double> pixes[6];
  vector<Ciphertext> pixes_en;

  Norm norm(encoder, encryptor, evaluator, relin_keys, parms);
  //encrype the originPic
  int K = parms.getk();

  //modify the size of y;
  if (y.size() % K != 0) {
    y.resize(y.size() - y.size() % K + K);
  }


  for (int i = 0; i < y.size(); i++) {
    pixes[i % K].push_back(y[i]);
  }


  cout << "the process is encoding orginal Pictures:" << endl;
  double dur;
  clock_t st, end;
  st = clock();
  Plaintext y_plain;
  Ciphertext y_en;
  for (int i = 0; i < K; i++) {
    encoder.encode(pixes[i], parms.getScale(), y_plain);
    encryptor.encrypt(y_plain,y_en);
    source.pushCipher(y_en);
  }

  pixes_en = source.getSecEn();
  Plaintext tem1;
  Ciphertext tem2;
  //encoder and encryted the index of X
  for (int i = 0; i < shares.size(); i++) {
    encoder.encode(shares[i].X[0], parms.getScale(), tem1);
    encryptor.encrypt(tem1, tem2);
    shares[i].x_en.push_back(tem2);
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.encode0 = (dur / CLOCKS_PER_SEC);
  cout << "the encoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;

  cout << "the Process-generate Shares by FHE is beginning:" << endl;
  st = clock();
  vector <Ciphertext>&temy = pixes_en;
  for (int pn = 0; pn < shares.size(); pn++) {
    Ciphertext Tem = shares[pn].generateCKKSShares(temy, evaluator, relin_keys, context,norm.ONES);
   }
  end = clock();
  dur = (double)(end - st);
  finalRes.generate0 = (dur / CLOCKS_PER_SEC);
  cout << "the Getting Shares' Process has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;


  cout << "the state of decoding the shares is beginning:" << endl;
  st = clock();
  //show the image
  for (int i = 0; i < shares.size(); i++) {
    if(printAns)
      cout << "the " << i << "th shares' pixes are:";
    shares[i].showShareCKKS(decryptor, encoder, printAns);
  }

  end = clock();
  dur = (double)(end - st);
  finalRes.decode0 = (dur / CLOCKS_PER_SEC);
  cout << "the state of decoding the shares has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
}


Picture recoryShareCKKS(Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& oriPic, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  cout << "--------------------------------------------\n---------------------------------------\n";
  cout << "Start Fully recovery the pic by FHE;\n";
  cout << "Gennerate the scheme base information:\n";

  EncryptionParameters params(scheme_type::CKKS);
  size_t poly_modulus_degree = degree;
  params.set_poly_modulus_degree(poly_modulus_degree);
  vector<int>coeff_modulus;

  coeff_modulus.push_back(min(picParms.getBaseLenOfModulus() + 20, 60));
  for (int i = 0; i < picParms.getMaxLevel(); i++) {
    coeff_modulus.push_back(picParms.getBaseLenOfModulus());
  }
  coeff_modulus.push_back(min(picParms.getBaseLenOfModulus() + 20, 60));

  params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus));

  auto context = SEALContext::Create(params);
  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  CKKSEncoder encoder(context);
  Plaintext tem;

  double dur;
  clock_t start, end;

  {
    cout<<encoder.slot_count();
  }

  {
    cout << "Start to encode the local shares:\n";
    cout << "waiting...." << endl;
    start = clock();

    for (int i = 0; i < uploadShares.size(); i++) {
      encoder.encode(uploadShares[i].X[0], picParms.getScale(), uploadShares[i].X_plain);
      encryptor.encrypt(uploadShares[i].X_plain, uploadShares[i].x_en[0]);
      uploadShares[i].temVec.clear();
      // using all the pixes to get a vector
      for (int j = 0; j < uploadShares[i].fx.size(); j++) {
        uploadShares[i].temVec.push_back((double)uploadShares[i].fx[j]);
      }
      encoder.encode(uploadShares[i].temVec, picParms.getScale(), uploadShares[i].fx_plain);
      encryptor.encrypt(uploadShares[i].fx_plain, uploadShares[i].fx_en[0]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "encoding procedure is finished.";
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
    finalRes.encode1 = (dur / CLOCKS_PER_SEC);
  }

  Norm norm(encoder,encryptor,evaluator,relin_keys,picParms);
  vector<Ciphertext>X, Ys;
  Plaintext invKT_plain;
  Ciphertext invKT_en;
  Picture recon_pic(picParms);
  
  {
    cout << "Start reconstructing Pics procedure:" << endl;
    start = clock();
    /// <summary>
    /// print invKT
    /// </summary>
    /// <param name="picParms"></param>
    /// <param name="uploadShares"></param>
    /// <param name="tools"></param>
    /// <param name="finalRes"></param>
    /// <param name="degree"></param>
    /// <param name="p_mo"></param>
    /// <returns></returns>
    //cout << "the recovery invKT is " << tools.invKT[0] << endl;
    encoder.encode(tools.invKT[0], picParms.getScale(), invKT_plain);
    encryptor.encrypt(invKT_plain, invKT_en);
    for (int i = 0; i < uploadShares.size(); i++) {
      X.push_back(uploadShares[i].x_en[0]);
      Ys.push_back(uploadShares[i].fx_en[0]);
    }
    vector<Ciphertext> a = fullRecoveryCKKS(X, Ys, invKT_en, encryptor, evaluator, encoder, decryptor, relin_keys, norm, picParms, context);
    for (int i = 0; i < a.size(); i++) {
      recon_pic.pushCipher(a[i]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "reconstructing Pics procedure is finished." << endl << endl;
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl;
    finalRes.generate1 = (dur / CLOCKS_PER_SEC);
  }

  recon_pic.setPic(recon_pic.DecryPicCKKS(decryptor, encoder, evaluator, norm, finalRes, printAns));
  recon_pic.compare(oriPic,norm);
  
  return recon_pic;
}

Picture recoryShareCKKS2(Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& oriPic, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  cout << "--------------------------------------------\n---------------------------------------\n";
  cout << "Start Fully recovery the pic by FHE;\n";
  cout << "Gennerate the scheme base information:\n";

  EncryptionParameters params(scheme_type::CKKS);
  size_t poly_modulus_degree = degree;
  params.set_poly_modulus_degree(poly_modulus_degree);
  vector<int>coeff_modulus;

  coeff_modulus.push_back(min(picParms.getBaseLenOfModulus() + 20, 60));
  for (int i = 0; i < picParms.getMaxLevel(); i++) {
    coeff_modulus.push_back(picParms.getBaseLenOfModulus());
  }
  coeff_modulus.push_back(min(picParms.getBaseLenOfModulus() + 20, 60));

  params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_modulus));

  auto context = SEALContext::Create(params);
  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  CKKSEncoder encoder(context);
  Plaintext tem;

  double dur;
  clock_t start, end;

  {
    cout << encoder.slot_count();
  }

  {
    cout << "Start to encode the local shares:\n";
    cout << "waiting...." << endl;
    start = clock();

    for (int i = 0; i < uploadShares.size(); i++) {
      uploadShares[i].temVec.clear();
      // using all the orignal to get a vector
      for (int j = 0; j < uploadShares[i].fx.size(); j++) {
        uploadShares[i].temVec.push_back((double)uploadShares[i].fx[j]+picParms.getP()*uploadShares[i].fx_div[j]);
      }
      encoder.encode(uploadShares[i].temVec, picParms.getScale(), uploadShares[i].fx_plain);
      encryptor.encrypt(uploadShares[i].fx_plain, uploadShares[i].fx_en[0]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "encoding procedure is finished.";
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
    finalRes.encode1 = (dur / CLOCKS_PER_SEC);
  }

  Norm norm(encoder, encryptor, evaluator, relin_keys, picParms);
  vector<double>X;
  vector<Ciphertext>Ys;
  Plaintext invKT_plain;
  Ciphertext invKT_en;
  Picture recon_pic(picParms);

  {
    cout << "Start reconstructing Pics procedure:" << endl;
    start = clock();

    encoder.encode(tools.invKT[0], picParms.getScale(), invKT_plain);
    encryptor.encrypt(invKT_plain, invKT_en);
    for (int i = 0; i < uploadShares.size(); i++) {
      X.push_back(uploadShares[i].X[0]);
      Ys.push_back(uploadShares[i].fx_en[0]);
    }
    vector<Ciphertext> a = fullRecoveryCKKS2(X, Ys, invKT_en, encryptor, evaluator, encoder, decryptor, relin_keys, norm, picParms, context);
    for (int i = 0; i < a.size(); i++) {
      recon_pic.pushCipher(a[i]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "reconstructing Pics procedure is finished." << endl << endl;
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl;
    finalRes.generate1 = (dur / CLOCKS_PER_SEC);
  }

  recon_pic.setPic(recon_pic.DecryPicCKKS(decryptor, encoder, evaluator, norm, finalRes, printAns));
  recon_pic.compare(oriPic, norm);

  return recon_pic;
}


void getShareByBFV(size_t x, Picture& source, vector<SharePic>& shares,
  Params& parms, Result& finalRes, bool printAns = true) {
  cout << "Start to generate secret shares by FHE(BFV):" << endl;
  vector<Ciphertext> X_en;
  EncryptionParameters params(scheme_type::BFV);
  size_t poly_modulus_degree = x;
  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 40));
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  auto context = SEALContext::Create(params);

  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  BatchEncoder encoder(context);
  size_t slot_count = encoder.slot_count();
  size_t row_size = slot_count / 2;

  vector<double> y = source.getSecFromDouble();
  vector<int64_t> pixes[6];
  vector<Ciphertext> pixes_en;

  Norm norm(encoder, encryptor, evaluator, relin_keys, parms);
  //encrype the originPic
  int K = parms.getk();

  //modify the size of y;
  if (y.size() % K != 0) {
    y.resize(y.size() - y.size() % K + K);
  }


  for (int i = 0; i < y.size(); i++) {
    pixes[i % K].push_back(y[i]);
  }
  int x_size = pixes[0].size();
  // normal the size of pixes
  for (int i = 0; i < 6; i++) {
    pixes[i].resize(slot_count);
  }

  cout << "the process is encoding orginal Pictures:" << endl;
  double dur;
  clock_t st, end;
  st = clock();
  Plaintext y_plain;
  Ciphertext y_en;
  for (int i = 0; i < K; i++) {
    encoder.encode(pixes[i], y_plain);
    encryptor.encrypt(y_plain, y_en);
    source.pushCipher(y_en);
  }

  pixes_en = source.getSecEn();
  Plaintext tem1;
  Ciphertext tem2;
  vector<int64_t>X(slot_count, 0LL);
  //encoder and encryted the index of X
  for (int i = 0; i < shares.size(); i++) {
    for (int j = 0; j < x_size; j++)
      X[j] = shares[i].X[0];
    encoder.encode(X, tem1);
    encryptor.encrypt(tem1, tem2);
    shares[i].x_en.push_back(tem2);
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.encode0 = (dur / CLOCKS_PER_SEC);
  cout << "the encoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;

  cout << "the Process-generate Shares by FHE is beginning:" << endl;
  st = clock();
  vector <Ciphertext>& temy = pixes_en;

  for (int pn = 0; pn < shares.size(); pn++) {
    Ciphertext Tem = shares[pn].generateBFVShares(temy, evaluator, relin_keys, context, norm.one_en);
  }
  end = clock();
  dur = (double)(end - st);
  finalRes.generate0 = (dur / CLOCKS_PER_SEC);
  cout << "the Getting Shares' Process has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;


  cout << "the state of decoding the shares is beginning:" << endl;
  st = clock();
  //show the image
  for (int i = 0; i < shares.size(); i++) {
    if (printAns)
      cout << "the " << i << "th shares' pixes are:";
    shares[i].showShareBFV(decryptor, encoder, printAns);
  }

  end = clock();
  dur = (double)(end - st);
  finalRes.decode0 = (dur / CLOCKS_PER_SEC);
  cout << "the state of decoding the shares has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
}

Picture recoryShareBFV(Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& oriPic, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  cout << "--------------------------------------------\n---------------------------------------\n";
  cout << "Start Fully recovery the pic by FHE;\n";
  cout << "Gennerate the scheme base information:\n";

  EncryptionParameters params(scheme_type::BFV);
  size_t poly_modulus_degree = degree;
  params.set_poly_modulus_degree(poly_modulus_degree);
  params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 40));
  params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
  auto context = SEALContext::Create(params);

  print_parameters(context);
  KeyGenerator kengen(context);
  PublicKey pkt = kengen.public_key();
  SecretKey skt = kengen.secret_key();
  RelinKeys relin_keys = kengen.relin_keys();
  Encryptor encryptor(context, pkt);
  Evaluator evaluator(context);
  Decryptor decryptor(context, skt);
  BatchEncoder encoder(context);
  Plaintext tem;
  size_t slot_count = encoder.slot_count();
  size_t row_size = slot_count / 2;

  double dur;
  clock_t start, end;

  {
    cout << encoder.slot_count() << endl;
  }

  {
    cout << "Start to encode the local shares:\n";
    cout << "waiting...." << endl;
    start = clock();

    // the number of Xid and responsed y pix
    vector<int64_t> Xtem(slot_count, 0LL), Ytem(slot_count, 0LL);
    for (int i = 0; i < uploadShares.size(); i++) {
      for (int j = 0; j < slot_count; j++) {
        Xtem[j] = uploadShares[i].X[0];
      }
      encoder.encode(Xtem, uploadShares[i].X_plain);
      encryptor.encrypt(uploadShares[i].X_plain, uploadShares[i].x_en[0]);
      uploadShares[i].temVec.clear();
      // using all the pixes to get a vector
      for (int j = 0; j < uploadShares[i].fx.size(); j++) {
        Ytem[j] = uploadShares[i].fx[j];
      }
      encoder.encode(Ytem, uploadShares[i].fx_plain);
      encryptor.encrypt(uploadShares[i].fx_plain, uploadShares[i].fx_en[0]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "encoding procedure is finished.";
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
    finalRes.encode1 = (dur / CLOCKS_PER_SEC);
  }

  Norm norm(encoder, encryptor, evaluator, relin_keys, picParms);
  vector<Ciphertext>X, Ys;
  Plaintext invKT_plain;
  Ciphertext invKT_en;
  Picture recon_pic(picParms);

  {
    cout << "Start reconstructing Pics procedure:" << endl;
    start = clock();
    /// <summary>
    /// print invKT
    /// </summary>
    /// <param name="picParms"></param>
    /// <param name="uploadShares"></param>
    /// <param name="tools"></param>
    /// <param name="finalRes"></param>
    /// <param name="degree"></param>
    /// <param name="p_mo"></param>
    /// <returns></returns>
    //cout << "the recovery invKT is " << tools.invKT[0] << endl;
    vector<int64_t> invK(slot_count, 0LL);
    for (int i = 0; i < slot_count; i++) {
      invK[i] = tools.invKT[0];
    }
    encoder.encode(invK, invKT_plain);
    encryptor.encrypt(invKT_plain, invKT_en);
    for (int i = 0; i < uploadShares.size(); i++) {
      X.push_back(uploadShares[i].x_en[0]);
      Ys.push_back(uploadShares[i].fx_en[0]);
    }
    vector<Ciphertext> a = fullRecoveryBFV(X, Ys, invKT_en, encryptor, evaluator, encoder, decryptor, relin_keys, norm, picParms, context);
    for (int i = 0; i < a.size(); i++) {
      recon_pic.pushCipher(a[i]);
    }
    end = clock();
    dur = (double)(end - start);
    cout << "reconstructing Pics procedure is finished." << endl << endl;
    cout << "Using time: " << (dur / CLOCKS_PER_SEC) << "s" << endl;
    finalRes.generate1 = (dur / CLOCKS_PER_SEC);
  }

  recon_pic.setPic(recon_pic.DecryPicBFV(decryptor, encoder, evaluator, norm, finalRes, printAns));
  recon_pic.compare(oriPic, norm);

  return recon_pic;
}
