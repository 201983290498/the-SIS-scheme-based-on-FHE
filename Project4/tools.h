#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>
#include<string>
#include<seal/seal.h>
using namespace std;
using namespace seal;

struct Error {
  double total;
  double mean;
  double var;
  double maxn;
  int dataScale;
};


/// <summary>
/// print vector
/// </summary>
/// <param name="vec"></param>
/// <param name="num"></param>
void print_vec(vector<double>&vec,int num,int pres=4) {
  ios old_fmt(nullptr);
  old_fmt.copyfmt(cout);
  cout << fixed << setprecision(pres);

  int len = num * 2;
  if (vec.size() <= len) {
    vec.resize(len + 1);
  }
  cout << "[";
  for (int i = 0; i < num; i++) {
    cout << vec[i] << ", ";
  }
  cout << "..., ";
  for (int i = vec.size() - num; i < vec.size(); i++) {
    cout << vec[i] << ", ";
  }
  cout << "]" << endl;
  cout.copyfmt(old_fmt);
}

void print_vec(vector<int>& vec, int num, int pres = 4) {
  ios old_fmt(nullptr);
  old_fmt.copyfmt(cout);
  cout << fixed << setprecision(pres);

  int len = num * 2;
  if (vec.size() <= len) {
    vec.resize(len + 1);
  }
  cout << "[";
  for (int i = 0; i < num; i++) {
    cout << vec[i] << ", ";
  }
  cout << "..., ";
  for (int i = vec.size() - num; i < vec.size(); i++) {
    cout << vec[i] << ", ";
  }
  cout << "]" << endl;
  cout.copyfmt(old_fmt);
}

void print_line(int tag = 1) {
  for (int i = 0; i < tag; i++) {
    cout << endl;
  }
}

void PrintDeVec(Decryptor& decryptor, CKKSEncoder& encoder, Ciphertext& en) {
  Plaintext plain;
  decryptor.decrypt(en, plain);
  vector<double> temAns;
  encoder.decode(plain,temAns);
  print_vec(temAns, 5, 4);
  print_line();
}

void PrintDeVec(Decryptor& decryptor, CKKSEncoder& encoder, vector<Ciphertext>& en, int tag = 1,int num=4) {
  Plaintext plain;
  vector<double> temAns;
  for (int i = 0; i < en.size(); i++) {
    decryptor.decrypt(en[i], plain);
    encoder.decode(plain, temAns);
    if (tag) {
      print_vec(temAns, num);
    }
    else {
      print_vec(temAns, 1);
    }
  }
  print_line();
}

Error cal_err(vector<double>& vec) {
  Error err;
  err.maxn = err.var = err.total = 0;
  for (int i = 0; i < vec.size(); i++) {
    err.total += vec[i];
    err.maxn = max(fabs((double)vec[i]), err.maxn);
  }
  err.mean = err.total / vec.size();
  for (int i = 0; i < vec.size(); i++) {
    err.var += pow(vec[i]-err.mean, 2);
  }
  err.var /= vec.size();
  err.dataScale = vec.size();
  return err;
}

Error cal_err(vector<int> vec) {
  Error err;
  err.maxn = err.var = err.total = 0;
  for (int i = 0; i < vec.size(); i++) {
    err.total += vec[i];
    err.maxn = max((double)vec[i], err.maxn);
  }
  err.mean = err.total / vec.size();
  for (int i = 0; i < vec.size(); i++) {
    err.var += pow(vec[i] - err.mean, 2);
  }
  err.var /= vec.size();
  err.dataScale = vec.size();
  return err;
}

void printError(Error& err) {
  cout << endl;
  ios old_fmt(nullptr);
  old_fmt.copyfmt(cout);
  cout << fixed << setprecision(3);
  cout << "the error detials[siz,totalsum, max, mean, var]: [" << err.dataScale << ", " << err.total << ", " << err.maxn << ", " << err.mean << ", " << err.var << "]" << endl;
  cout.copyfmt(old_fmt);
}

void outPutCiphertext(vector<Ciphertext>& text, string filePre) {

  cout << "the coeff_mod_count is: " << text[0].coeff_mod_count() <<
    "\n the poly_modulus_degree is: " << text[0].poly_modulus_degree()<<endl;

  for (int i = 0; i < text.size(); i++) {
    string temName = filePre + to_string(i) + ".txt";
    fstream out(temName, ios::out|ios::binary);
    text[i].save(out);
    out.close();
  }
}

void outPutCiphertext(vector<Plaintext>& text, string filePre) {
  cout << filePre << endl;
  for (int i = 0; i < text.size(); i++) {
    cout << "the coeff_count is: " << text[0].coeff_count() <<
      "\n the poly_modulus_degree is: " << text[0].nonzero_coeff_count() << endl;
    string temName = filePre + to_string(i) + ".txt";
    fstream out(temName, ios::out | ios::binary);
    text[i].save(out);
    out.close();
  }
}
