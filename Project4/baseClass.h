#pragma once
#include"tools.h"
#include<fstream>
#include"mymath.h"
#include"math_test.h"

#define ll long long
using namespace std;
using namespace seal;


class Result {
public:
  Result() {

  }
  void writeResult(string filename="result.txt") {
    fstream file(filename, ios::out | ios::app);
    file << l << "," << w << "," << t << "," << n << "," << p << "," << x_range << ",";
    file << poly_d0 << "," << plain_m0 << ",";
    file << poly_d1 << "," << plain_m1 << ",";
    file << encode0 << "," << generate0 << "," << decode0 << ",";
    file << encode1 << "," << generate1 << "," << decode1 << "," << correct << endl;
    file.close();
  }
  double encode0, generate0, decode0;
  double encode1, generate1, decode1;
  double correct;
  int l, w, t, n, p, x_range;
  int poly_d0, plain_m0, poly_d1, plain_m1;
};

class Params {

public:
  friend class Picture;
  friend class SharePic;
  friend class DecTools;
  /**
    the constuctor
  */
  Params() {

  }

  Params(int l = 2, int w = 2, int p = 251, int n = 6, int k = 4, int x = 10){
    L = l;
    W = w;
    ModP = p;
    this->n = n;
    this->k = k;
    rangeX = x;
    baseLenOfModulus = 40;
    this->batches = ((L * W % k == 0) ? L * W / k : (L * W - L * W % k + k) / k);
  }
  int* getSize() {
    int siz[] = { L,W };
    return siz;
  }

  int getP() {
    return ModP;
  }
  void setP(int  P) {
    ModP = P;
  }

  int getN() {
    return n;
  }
  void setN(int n){
    this->n = n;
  }

  int getk() {
    return k;
  }
  void setK(int K) {
    k = K;
  }

  int getXrange() {
    return rangeX;
  }
  void setXrange(int range) {
    rangeX = range;
  }

  void setMaxLevel(int x) {
    maxLevel = x;
  }
  int getMaxLevel() {
    return maxLevel;
  }

  void setScale(double scale) {
    this->scale = scale;
  }
  double getScale() {
    return this->scale;
  }

  int getBaseLenOfModulus() {
    return baseLenOfModulus;
  }
  void setBaseLenOfModulus(int len) {
    baseLenOfModulus = len;
  }


  int getBatches() {
    return batches;
  }
  void setBatch(int batches) {
    this->batches = batches;
  }

private:
  int L, W;
  int ModP;
  int n;
  int k,batches;
  int rangeX;
  int maxLevel;
  double scale;
  int baseLenOfModulus;
};

class Norm {
public:
  Norm() {}
 //the scheme of BFV
  Norm(IntegerEncoder& encoder, Encryptor& encryptor, Params& parms) {
    Plaintext zeros = encoder.encode(0), one = encoder.encode(1), P_splain = encoder.encode(parms.getP());
    encryptor.encrypt(zeros, zeros_en);
    encryptor.encrypt(one, one_en);
    encryptor.encrypt(P_splain, P_sen);
    modP = parms.getP();
  }

  //the scheme of CKKS
  Norm(CKKSEncoder& encoder, Encryptor& encryptor, Evaluator& evaluator,RelinKeys& relin_keys, Params& parms) {
    //generate the vector of ones
    Ciphertext tem;
    modP = parms.getP();
    for (int i = 0; i <= parms.getMaxLevel(); i++) {
      ONES.push_back(tem);
    }
    Plaintext one, zero, p_plain;
    encoder.encode(1.0,parms.getScale(),one);
    encoder.encode(0,parms.getScale(),zero);
    encoder.encode(parms.getP(),parms.getScale(),p_plain);
    encryptor.encrypt(one, tem);
    encryptor.encrypt(zero, ZERO);
    encryptor.encrypt(p_plain, P_sen);
    for (int i = parms.getMaxLevel(); i >= 0; i--) {
      ONES[i] = tem;
      if (i != 0) {
        evaluator.multiply_inplace(tem, tem);
        evaluator.relinearize_inplace(tem, relin_keys);
        evaluator.rescale_to_next_inplace(tem);
      }
    }
  }

  // the scheme of BFVBatchProcess
  Norm(BatchEncoder& encoder, Encryptor& encryptor, Evaluator& evaluator, RelinKeys& relin_keys, Params& parms) {
    //generate the vector of ones
    // the vector is full of
    //the same number 0, 1, modP
    Ciphertext tem;
    modP = parms.getP();
    size_t slot_count = encoder.slot_count();
    size_t row_size = slot_count / 2;

    vector<int64_t> ONE(slot_count, 0LL);
    vector<int64_t> ZERO(slot_count, 0LL);
    vector<int64_t> P(slot_count, 0LL);
    for (int i = 0; i < slot_count; i++) {
      ONE[i] = 1LL;
      P[i] = modP;
    }
    Plaintext one, zero, p_plain;
    encoder.encode(ONE, one);
    encoder.encode(ZERO, zero);
    encoder.encode(P, p_plain);
    encryptor.encrypt(one, one_en);
    encryptor.encrypt(zero, this->ZERO);
    encryptor.encrypt(p_plain, P_sen);
  }


  Ciphertext zeros_en, P_sen, one_en;
  vector<Ciphertext> ONES;
  Ciphertext ZERO;
  int modP;

  // the error the ckks
  Error virError, realError;
  vector<double> virErrorData, realErrorData;
  vector<int> realErrorIndex;

};

class Picture {
  friend class SharePic;
public:
  Picture(){}
  Picture(Params& parms) {
    siz[0] = parms.L;
    siz[1] = parms.W;
    for (int i = 0; i < siz[0] * siz[1]; i++) {
      origin_pic.push_back(0);
      origin_sec.push_back(0);
      aerfa.push_back(0);
    }
    batches = parms.getBatches();
  }

  vector<ll>& getSec() {
    return origin_sec;
  }

  vector<double> getSecFromDouble() {
    if (origin_secCKKS.size() == 0) {
      //change the type of orign_sec
      for (int i = 0; i < origin_sec.size(); i++) {
        origin_secCKKS.push_back((double)origin_sec[i]);
      }
    }
    return origin_secCKKS;
  }

  vector<Ciphertext> getSecEn() {
    return pic_en;
  }

  int* getSiz() {
    return siz;
  }

  /**
  1. get the originPic
  2. the default randMul of pix is 0
  3. we can control the input pics
  */
  void generatePic(Params& parms, int randRange = 1) {
    for (int i = 0; i < siz[0]; i++) {
      for (int j = 0; j < siz[1]; j++) {
        int index = i * parms.W + j;
        origin_pic[index] = rand() % parms.ModP;
        aerfa[index] = (ll)rand() % randRange;
        origin_sec[index] = aerfa[index] * parms.ModP + origin_pic[index];
      }
    }
  }

  void PrintParms(vector<ll>& vec) {
    for (int i = 0; i < vec.size(); i++)
      cout << vec[i] << " ";
    cout << endl;
  }

  void pushCipher(Ciphertext& text) {
    pic_en.push_back(text);
  }

  void printPic() {
    for (int i = 0; i < origin_sec.size(); i++) {
      cout << origin_pic[i] << " ";
    }
    cout << endl << endl;
  }

  void pushPies(ll pixes, int index) {
    origin_pic[index] = pixes;
  }

  vector<int> DecryPic(Decryptor& decryptor, IntegerEncoder& encoder, Evaluator& evaluator, Norm& norm, Result& finalRes, bool printAns = true) {
    Plaintext text1;
    vector<Ciphertext>a = pic_en;
    Ciphertext zeros_en = norm.zeros_en, en1;
    cout << "download from  the center server;\n";
    cout << "\n\n\nstart to decoder locally\n";
    double dur;
    clock_t st, end;
    st = clock();
    vector<int>ans_fl;
    BigUInt tem_de;
    BigUInt temP(100, norm.modP);
    BigUInt temAns;
    // cout << "the orgin data of tem_de:\n";
    for (int i = 0; i < a.size(); i++) {
      decryptor.decrypt(a[i], text1);
      try {
        tem_de = encoder.decode_biguint(text1);
        tem_de.divrem(temP, temAns);
        //cout << tem_de.to_dec_string() << endl;
        ans_fl.push_back(temAns.to_double());
      }
      catch (invalid_argument e) {
        evaluator.sub(zeros_en, a[i], en1);
        check1(decryptor, en1);
        decryptor.decrypt(en1, text1);
        try {
          tem_de = encoder.decode_biguint(text1);
          cout << tem_de.to_dec_string() << endl;
          tem_de.divrem(temP, temAns);
          temAns = temP - temAns;
          temAns.divrem(temP, temAns);
          ans_fl.push_back(temAns.to_double());
        }
        catch (invalid_argument e) {
          ans_fl.push_back(0);
        }
      }

    }
    end = clock();
    dur = (double)(end - st);
    finalRes.encode0 = (dur / CLOCKS_PER_SEC);
    cout << "the decoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;

    if (printAns) {
      cout << "the recovery result of pixes using FHE:\n";
      for (int i = 0; i < ans_fl.size(); i++) {
        cout << ans_fl[i] << " ";
      }
      cout << endl;
    }
    return ans_fl;
  }

  vector<int> DecryPicCKKS(Decryptor& decryptor, CKKSEncoder& encoder, Evaluator& evaluator, Norm& norm, Result& finalRes, bool printAns = true) {
    Plaintext text1;
    vector<Ciphertext>a = pic_en;
    cout << a.size() << endl;
    vector<double> temAns[10];
    vector<int>ans_fl;
    Ciphertext zeros_en = norm.zeros_en, en1;
    cout << "download from  the center server;\n";
    cout << "\n\n\nstart to decoder locally\n";
    {
      double dur;
      clock_t st, end;
      st = clock();

      // cout << "the orgin data of tem_de:\n";
      for (int i = 0; i < a.size(); i++) {
        decryptor.decrypt(a[i], text1);
        encoder.decode(text1, temAns[i]);
      }
      origin_pic.clear();
      cout << norm.modP << endl;
      for (int i = 0; i < batches; i++)
        for (int j = 0; j < a.size(); j++) {
          ans_fl.push_back((around(temAns[j][i]) % norm.modP + norm.modP) % norm.modP);
          //ans_fl.push_back((around(temAns[j][i])));
          norm.virErrorData.push_back(around(temAns[j][i]) - temAns[j][i]);
        }
      end = clock();
      dur = (double)(end - st);
      finalRes.encode0 = (dur / CLOCKS_PER_SEC);
      cout << "the decoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
    }

    {
      norm.virError = cal_err(norm.virErrorData);
      printError(norm.virError);
    }

    if (printAns) {
      cout << "the recovery result of pixes using FHE:\n";
      for (int i = 0; i < ans_fl.size(); i++) {
        cout << ans_fl[i] << " ";
      }
      cout << endl;
    }
    return ans_fl;
  }

  vector<int> DecryPicBFV(Decryptor& decryptor, BatchEncoder& encoder, Evaluator& evaluator, Norm& norm, Result& finalRes, bool printAns = true) {
    Plaintext text1;
    vector<Ciphertext>a = pic_en;
    cout << a.size() << endl;
    vector<int64_t> temAns[10];
    vector<int>ans_fl;
    Ciphertext zeros_en = norm.zeros_en, en1;
    cout << "download from  the center server;\n";
    cout << "\n\n\nstart to decoder locally\n";
    {
      double dur;
      clock_t st, end;
      st = clock();

      // cout << "the orgin data of tem_de:\n";
      for (int i = 0; i < a.size(); i++) {
        decryptor.decrypt(a[i], text1);
        encoder.decode(text1, temAns[i]);
      }
      origin_pic.clear();
      cout << norm.modP << endl;
      for (int i = 0; i < batches; i++)
        for (int j = 0; j < a.size(); j++) {
          ans_fl.push_back((around(temAns[j][i]) % norm.modP + norm.modP) % norm.modP);
          //ans_fl.push_back((around(temAns[j][i])));
          norm.virErrorData.push_back(around(temAns[j][i]) - temAns[j][i]);
        }
      end = clock();
      dur = (double)(end - st);
      finalRes.encode0 = (dur / CLOCKS_PER_SEC);
      cout << "the decoding processing has finished. The total time cost is:" << (dur / CLOCKS_PER_SEC) << "s" << endl << endl;
    }

    {
      norm.virError = cal_err(norm.virErrorData);
      printError(norm.virError);
    }


    if (printAns) {
      cout << "the recovery result of pixes using FHE:\n";
      for (int i = 0; i < ans_fl.size(); i++) {
        cout << ans_fl[i] << " ";
      }
      cout << endl;
    }

    return ans_fl;
  }

  void setBatches(int batches) {
    this->batches = batches;
  }

  int getBatches() {
    return this->batches;
  }

  void setPic(vector<int> a) {
    origin_pic = a;
  }

  double compare(const Picture& tem)const {
    double flag = 0, total = 0;
    if (tem.origin_pic.size() != origin_pic.size())
      flag = 0;
    for (int i = 0; i < origin_pic.size(); i++)
      if (tem.origin_pic[i] != origin_pic[i]) {
        total = total + 1;
      }
      else {
        total = total + 1;
        flag++;
      }
    cout << "the correct rate of ans is " << (flag ? (flag / total * 100) : 0) << "%" << endl;
    return (flag ? (flag / total * 100) : 0);
  }

  double compare(const Picture& tem,Norm& norm)const {

    if (tem.origin_pic.size() != origin_pic.size()) {
      cout << "ERROR!!";
      throw "sizes are different";
    }
    for (int i = 0; i < origin_pic.size(); i++)
      if (tem.origin_pic[i] != origin_pic[i]) {
        norm.realErrorData.push_back(fabs(tem.origin_pic[i] - origin_pic[i]));
        norm.realErrorIndex.push_back(i);
      }
    norm.realError = cal_err(norm.realErrorData);
    cout << endl;
    printError(norm.realError);
  }

private:
  vector<int> origin_pic;
  vector<ll> origin_sec, aerfa;
  vector<Ciphertext>pic_en;
  int siz[2] = { 0,0 };

  vector<double>origin_secCKKS;
  int batches;
};

//shares
class SharePic {

public:
  friend class Picture;
  SharePic(Params& parms) {
    range_x = parms.rangeX;
    P = parms.ModP;
  }

  // using Plain to get sharePixes
  void addNewPixByPlain(vector<ll>y, int x,bool BFV=true) {
    ll f = 0;
    int k = y.size();
    for (int i = 0; i < k; i++) {
      f += (long long)pow((ll)x, i) * y[i];
    }
    if (BFV) 
      X.push_back(x);
    fx_div.push_back((f - f % P) / P);
    fx.push_back((int)(f % P));
  }

  vector<ll> getX() {
    return X;
  }

  vector<int> getfx() {
    return fx;
  }

  void PrintFx() {
    for (int i = 0; i < fx.size(); i++) {
      cout << fx[i] << " ";
    }
    cout << endl << endl;
  }

  Ciphertext& addNewPixByCipher(vector<Ciphertext>y, Evaluator& evaluator, RelinKeys& rk, int index = 0) {
    Ciphertext tem = y[0], tem2;
    for (int i = 1; i < y.size(); i++) {
      evaluator.exponentiate(x_en[index], i, rk, tem2);
      evaluator.multiply_inplace(tem2, y[i]);
      evaluator.add_inplace(tem, tem2);
    }
    evaluator.relinearize_inplace(tem2, rk);
    fx_en.push_back(tem);
    return tem;
  }

  Ciphertext generateCKKSShares(vector<Ciphertext>y, Evaluator& evaluator, RelinKeys& rk, shared_ptr<seal::SEALContext>& context,vector<Ciphertext>& ONES) {
    int k = y.size();
    Ciphertext ans = y[0],tem1;
    for (int i = 1; i < k; i++) {
      vector<Ciphertext> tem;
      tem.push_back(y[i]);
      for (int j = 1; j <= i; j++)
        tem.push_back(x_en[0]);
      tem1 = cal_mut(tem, evaluator, rk, context,ONES);
      add(evaluator, ans, tem1, context);
    }
    fx_en.push_back(ans);
    return ans;
  }

  Ciphertext generateBFVShares(vector<Ciphertext>y, Evaluator& evaluator, RelinKeys& rk, shared_ptr<seal::SEALContext>& context, Ciphertext& ONES) {
    int k = y.size();
    Ciphertext ans = y[0], tem1;
    for (int i = 1; i < k; i++) {
      vector<Ciphertext> tem;
      tem.push_back(y[i]);
      for (int j = 1; j <= i; j++)
        tem.push_back(x_en[0]);
      tem1 = cal_mutBFV(tem, evaluator, rk, ONES);
      evaluator.add_inplace(ans, tem1);
    }
    fx_en.push_back(ans);
    return ans;
  }

  void showShare(Decryptor& decryptor, IntegerEncoder& encoder, bool printAns = true) {
    Plaintext tem;
    bool flag = true;
    for (int i = 0; i < fx_en.size(); i++) {
      decryptor.decrypt(fx_en[i], tem);
      int temnum = encoder.decode_int64(tem) % P;
      if(printAns)
        cout << temnum << " ";
      if (temnum != fx[i])
        flag = false;
    }
    if (flag)
      cout << "the ans is correct";
    else
      cout << "something is error!";
    cout << endl;
  }

  void showShareCKKS(Decryptor& decryptor, CKKSEncoder& encoder,bool printAns = true) {
    Plaintext tem;
    bool flag = true;
    vector<double> share_de;
    decryptor.decrypt(fx_en[0], tem);
    encoder.decode(tem, share_de);
    for (int i = 0; i < X_len; i++) {
      //cout << fixed << setprecision(2) << share_de[i] << " " << fabs(around(share_de[i]) - share_de[i]) << " ";
      if(printAns)
        cout << around(share_de[i]) % P << " ";
      if ( around(share_de[i]) % P != fx[i])
        flag = false;
    }
    if (flag)
      cout << "the ans is correct";
    else
      cout << "something is error!";
    cout << endl;
  }

  void showShareBFV(Decryptor& decryptor, BatchEncoder& encoder, bool printAns = true) {
    Plaintext tem;
    bool flag = true;
    vector<int64_t> share_de;
    decryptor.decrypt(fx_en[0], tem);
    encoder.decode(tem, share_de);
    for (int i = 0; i < X_len; i++) {
      //cout << fixed << setprecision(2) << share_de[i] << " " << fabs(around(share_de[i]) - share_de[i]) << " ";
      if (printAns)
        cout << around(share_de[i]) % P << " ";
      if (around(share_de[i]) % P != fx[i])
        flag = false;
    }
    if (flag)
      cout << "the ans is correct";
    else
      cout << "something is error!";
    cout << endl;
  }


  vector<int> fx,fx_div;
  vector<ll> X;
  int range_x;
  vector<Ciphertext> x_en;
  vector<Ciphertext> fx_en;
  ll pix_de;
  ll pix_sr;
  int rad;//the use of rad
  int P;

  // in the scheme CKKS,the X's number is 1;
  // however the hidden information is X_len;
  int X_len;

  Plaintext X_plain, fx_plain;
  vector<double>temVec;
};

//decrypted
class DecTools {
public:
  vector<long long>m[64*64];
  vector<long long>K[64 * 64];
  vector<long long>invm[64 * 64];
  vector<long long>A[64 * 64];
  long long KT[64 * 64];
  long long invKT[64 * 64];
};

