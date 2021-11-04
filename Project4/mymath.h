#pragma once
#include<cmath>
#include<time.h>
#include<vector>
#include "seal/seal.h"
#include"math_test.h"

using namespace std;
using namespace seal;

void check1(Decryptor& decryptor, Ciphertext xx) {
  int bits = decryptor.invariant_noise_budget(xx);
  if (bits > 0) {
    return;
  }
  else {
    cout << "the noise is out of range;\n";
    exit(0);
  }
}


long long mypow(long long base, int mi, int MOD) {
  long long ans = 1;
  while (mi) {
    if (mi % 2) {
      ans = ans * base % MOD;
    }
    mi = mi >> 1;
    base = base * base % MOD;
  }
  return ans;
}

vector<int> getD(vector<int>& X, int k) {
  vector<int> D[10];
  vector<int>ans;
  for (int i = 0; i < X.size(); i++)
    X[i] = -X[i];

  for (int i = 0; i < k; i++) {//from x0 to xk-1
    vector<int>b;

    for (int j = 0; j < k; j++) {//the range of generate b are b 0~k-2
      if (j == i)continue;
      b.push_back(X[j]);
    }
    long long presum[10] = { 0,0,0,0,0,0,0,0,0,0 };
    long long nsum[10] = { 0,0,0,0,0,0,0,0,0,0 };
    for (int j = 0; j < k; j++) {
      if (j == 0) {
        D[i].push_back(1);
      }
      else if (j == 1) {
        for (int ii = 1; ii < k; ii++) {
          presum[ii] = presum[ii - 1] + b[ii - 1];//0~k-2
        }
        D[i].push_back(presum[k - j]);//1~k-1
      }
      else {
        for (int ii = 1; ii <= k - j; ii++) {
          long long tem = b[ii - 1] * (presum[k - j + 1] - presum[ii]);//get the presum
          nsum[ii] = tem + nsum[ii - 1];
        }
        D[i].push_back(nsum[k - j]);
        swap(presum, nsum);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      ans.push_back(D[i][j]);
    }
  }
  return ans;
}

vector<double> getD(vector<double>& X, int k) {
  vector<int> D[10];
  vector<double>ans;
  for (int i = 0; i < X.size(); i++)
    X[i] = -X[i];

  for (int i = 0; i < k; i++) {//from x0 to xk-1
    vector<int>b;

    for (int j = 0; j < k; j++) {//the range of generate b are b 0~k-2
      if (j == i)continue;
      b.push_back(X[j]);
    }
    long long presum[10] = { 0,0,0,0,0,0,0,0,0,0 };
    long long nsum[10] = { 0,0,0,0,0,0,0,0,0,0 };
    for (int j = 0; j < k; j++) {
      if (j == 0) {
        D[i].push_back(1);
      }
      else if (j == 1) {
        for (int ii = 1; ii < k; ii++) {
          presum[ii] = presum[ii - 1] + b[ii - 1];//0~k-2
        }
        D[i].push_back(presum[k - j]);//1~k-1
      }
      else {
        for (int ii = 1; ii <= k - j; ii++) {
          long long tem = b[ii - 1] * (presum[k - j + 1] - presum[ii]);//get the presum
          nsum[ii] = tem + nsum[ii - 1];
        }
        D[i].push_back(nsum[k - j]);
        swap(presum, nsum);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      ans.push_back(D[i][j]);
    }
  }
  return ans;
}


vector<Ciphertext> getD(vector<Ciphertext>& X, Evaluator& evaluator, Decryptor& decryptor, RelinKeys relinkeys, int k, Ciphertext& ZERO, Ciphertext& ONE) {
  vector<Ciphertext> D[10];
  vector<Ciphertext>ans;

  for (int i = 0; i < X.size(); i++) {
    evaluator.sub(ZERO, X[i], X[i]);
  }

  for (int i = 0; i < k; i++) {//from x0 to xk-1
    vector<Ciphertext>b;

    for (int j = 0; j < k; j++) {//the range of generate b are b 0~k-2
      if (j == i)continue;
      b.push_back(X[j]);
    }
    Ciphertext presum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    Ciphertext nsum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    for (int j = 0; j < k; j++) {
      if (j == 0) {
        D[i].push_back(ONE);
      }
      else if (j == 1) {
        for (int ii = 1; ii < k; ii++) {
          evaluator.add(presum[ii - 1], b[ii - 1], presum[ii]);
        }
        evaluator.relinearize_inplace(presum[k - j], relinkeys);
        //check1(decryptor, presum[k - j]);
        D[i].push_back(presum[k - j]);//1~k-1
      }
      else {
        for (int ii = 1; ii <= k - j; ii++) {
          Ciphertext tem;
          evaluator.sub(presum[k - j + 1], presum[ii], tem);
          evaluator.multiply_inplace(tem, b[ii - 1]);
          evaluator.relinearize_inplace(tem, relinkeys);
          evaluator.add(tem, nsum[ii - 1], nsum[ii]);
        }
        evaluator.relinearize_inplace(nsum[k - j], relinkeys);
        //check1(decryptor, nsum[k - j]);
        D[i].push_back(nsum[k - j]);
        swap(presum, nsum);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      ans.push_back(D[i][j]);
    }
  }
  return ans;
}

vector<Ciphertext> getDCKKS(vector<Ciphertext>& X, Evaluator& evaluator, Decryptor& decryptor, RelinKeys relinkeys, int k, Ciphertext& ZERO, vector<Ciphertext>& ONES, shared_ptr<seal::SEALContext>& context) {
  vector<Ciphertext> D[10];
  vector<Ciphertext>ans;

  for (int i = 0; i < X.size(); i++) {
    evaluator.sub(ZERO, X[i], X[i]);
  }

  for (int i = 0; i < k; i++) {//from x0 to xk-1
    vector<Ciphertext>b;

    for (int j = 0; j < k; j++) {//the range of generate b are b 0~k-2
      if (j == i)continue;
      b.push_back(X[j]);
    }
    Ciphertext presum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    Ciphertext nsum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    for (int j = 0; j < k; j++) {
      if (j == 0) {
        D[i].push_back(ONES[ONES.size()-1]);
      }
      else if (j == 1) {
        for (int ii = 1; ii < k; ii++) {
          add(evaluator, presum[ii - 1], b[ii - 1], presum[ii], context);
        }
        //check1(decryptor, presum[k - j]);
        D[i].push_back(presum[k - j]);//1~k-1
      }
      else {
        for (int ii = 1; ii <= k - j; ii++) {
          Ciphertext tem;
          sub(evaluator,presum[k - j + 1], presum[ii], tem,context);
          tem = cal_mut({ tem, b[ii - 1] },evaluator,relinkeys,context,ONES);
          add(evaluator,tem, nsum[ii - 1], nsum[ii],context);
        }
        D[i].push_back(nsum[k - j]);
        swap(presum, nsum);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      ans.push_back(D[i][j]);
    }
  }
  return ans;
}

vector<Ciphertext> getDBFV(vector<Ciphertext>& X, Evaluator& evaluator, Decryptor& decryptor, RelinKeys relinkeys, int k, Ciphertext& ZERO, Ciphertext ONES) {
  vector<Ciphertext> D[10];
  vector<Ciphertext>ans;

  for (int i = 0; i < X.size(); i++) {
    evaluator.sub(ZERO, X[i], X[i]);
  }

  for (int i = 0; i < k; i++) {//from x0 to xk-1
    vector<Ciphertext>b;

    for (int j = 0; j < k; j++) {//the range of generate b are b 0~k-2
      if (j == i)continue;
      b.push_back(X[j]);
    }
    Ciphertext presum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    Ciphertext nsum[10] = { ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO };
    for (int j = 0; j < k; j++) {
      if (j == 0) {
        D[i].push_back(ONES);
      }
      else if (j == 1) {
        for (int ii = 1; ii < k; ii++) {
          evaluator.add(presum[ii - 1], b[ii - 1], presum[ii]);
        }
        //check1(decryptor, presum[k - j]);
        D[i].push_back(presum[k - j]);//1~k-1
      }
      else {
        for (int ii = 1; ii <= k - j; ii++) {
          Ciphertext tem;
          evaluator.sub(presum[k - j + 1], presum[ii], tem);
          evaluator.multiply_inplace(tem, b[ii - 1]);
          evaluator.relinearize_inplace(tem, relinkeys);
          evaluator.add(tem, nsum[ii - 1], nsum[ii]);
        }
        D[i].push_back(nsum[k - j]);
        swap(presum, nsum);
      }
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      ans.push_back(D[i][j]);
    }
  }
  return ans;
}



Ciphertext mypow(Evaluator& evaluator, Ciphertext baseTE, Decryptor& decryptor, int mi, Ciphertext ONE, RelinKeys relinkeys, IntegerEncoder& encoder) {//the important shortcut is the number is too big to be a effective scheme;
  Ciphertext ans = ONE;
  Plaintext plain;
  cout << "trying to cal the  invK\n";
  while (mi) {
    cout << "Noise budget of invK:" << decryptor.invariant_noise_budget(baseTE) << "bits" << endl;
    if (mi % 2) {
      evaluator.multiply_inplace(ans, baseTE);
      evaluator.relinearize_inplace(ans, relinkeys);
    }
    mi = mi >> 1;
    evaluator.multiply_inplace(baseTE, baseTE);
    evaluator.relinearize_inplace(baseTE, relinkeys);
  }
  return ans;
}
