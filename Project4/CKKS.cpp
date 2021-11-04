#include<string>  
#include "math_test.h"
#include"PlainProcess.h"
#include"cipherProcess.h"
#include<map>

using namespace std;
using namespace seal;

void functionTestBFV(int mode, Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& pic,int writeSimpleData, string filename, ifstream& datain, vector<SharePic>& shares,int degree = 8192, int p_mo = 1024, bool printAns = true) {
  int outPutCip;
  cout << "whether to output the CiphereText(0|1):";
  cin >> outPutCip;
  if (mode == 1) {
    Picture recon_pic = recoryShareBFV(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
    finalRes.correct = recon_pic.compare(pic);
    cout << "the invKT is: " << tools.invKT[0] << endl;
    if (writeSimpleData) {
      finalRes.writeResult(filename);
    }
  }
  else if (mode == 2) {
    cout << "please input the range of coeff_modulus bits[start, step, end]";
    int start, step, end;
    cin >> start >> step >> end;
    for (int i = start; i <= end; i += step) {
      Picture recon_pic = recoryShareBFV(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
      finalRes.correct = recon_pic.compare(pic);
      cout << "the invKT is: " << tools.invKT[0] << endl;
      if (writeSimpleData) {
        finalRes.writeResult(filename);
      }
    }
  }
  else if (mode == 3) {
    int ti;
    cout << "please input the repeat times:";
    cin >> ti;
    while (ti--) {
      Picture recon_pic = recoryShareBFV(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
      finalRes.correct = recon_pic.compare(pic);
      cout << "the invKT is: " << tools.invKT[0] << endl;
      if (writeSimpleData) {
        finalRes.writeResult(filename);
      }
    }
  }
  else if (mode == 4) {
    map<vector<int>, double> cal_ans;
    int bits;
    cout << "please input the repeat times and coeff_modulus bits:";
    cin >> bits;
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        vector<int> recovery_ids;
        vector<SharePic>UploadShares;
        for (int k = 0; k < 6; k++)
          if (k != i && k != j) {
            recovery_ids.push_back(k);
            UploadShares.push_back(shares[k]);
          }
        recoveryByPlainCKKS(UploadShares, picParms, tools);
        if (tools.invKT[0] > picParms.getP() / 2) {
          tools.invKT[0] -= picParms.getP();
        }
        Picture recon_pic = recoryShareBFV(picParms, UploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
        finalRes.correct = recon_pic.compare(pic);
        cout << "the invKT is: " << tools.invKT[0] << endl;
        cal_ans[recovery_ids] = finalRes.correct;
      }
    }
    fstream out("all.txt", ios::out);
    map<vector<int>, double>::iterator iter = cal_ans.begin();
    while (iter != cal_ans.end()) {
      vector<int> tem = iter->first;
      for (int i = 0; i < tem.size(); i++) {
        out << tem[i] << ",";
      }
      out << iter->second << endl;
      iter++;
    }
    out.close();
  }

  if (outPutCip)
  {
    Plaintext tem;
    vector<Ciphertext>x, fx;
    vector<Plaintext>x_plain, fx_plain;
    for (int i = 0; i < uploadShares.size(); i++) {
      x.push_back(uploadShares[i].x_en[0]);
      x_plain.push_back(uploadShares[i].X_plain);
      fx.push_back(uploadShares[i].fx_en[0]);
      fx_plain.push_back(uploadShares[i].fx_plain);
    }
    tem = x_plain[0];
    for (int i = 0; i < 10; i++) {
      cout << tem[i] << " ";
    }
    tem = fx_plain[0];
    
    for (int i = 0; i < 10; i++) {
      cout << tem[i] << " ";
    }
    cout << endl;
    cout << "output the shared pic:\n";
    outPutCiphertext(fx, "ciphertext");
    cout << "output the x_id:\n";
    outPutCiphertext(x, "xid");
    cout << "output the plainText of x\n";
    outPutCiphertext(x_plain, "x_plaintext");
    cout << "output the plainText of vector y";
    outPutCiphertext(fx_plain, "vector_fx");
  }
}

void functionTest(int mode, Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& pic, int writeSimpleData, string filename, ifstream& datain, vector<SharePic>& shares, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  int outPutCip;
  cout << "whether to output the CiphereText(0|1):";
  cin >> outPutCip;
  if (mode == 1) {
    //每组x,y输入，然后就可以得到想要的y
//recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
    picParms.setBaseLenOfModulus(50);
    picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
    Picture recon_pic = recoryShareCKKS(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
    finalRes.correct = recon_pic.compare(pic);
    cout << "the invKT is: " << tools.invKT[0] << endl;
    if (writeSimpleData) {
      finalRes.writeResult(filename);
    }
  }
  else if (mode == 2) {
    cout << "please input the range of coeff_modulus bits[start, step, end]";
    int start, step, end;
    cin >> start >> step >> end;
    for (int i = start; i <= end; i += step) {
      picParms.setBaseLenOfModulus(i);
      picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
      Picture recon_pic = recoryShareCKKS(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
      finalRes.correct = recon_pic.compare(pic);
      cout << "the invKT is: " << tools.invKT[0] << endl;
      if (writeSimpleData) {
        finalRes.writeResult(filename);
      }
    }
  }
  else if (mode == 3) {
    int ti, bits;
    cout << "please input the repeat times and coeff_modulus bits:";
    cin >> ti >> bits;
    while (ti--) {
      //每组x,y输入，然后就可以得到想要的y
      //recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
      picParms.setBaseLenOfModulus(bits);
      picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
      Picture recon_pic = recoryShareCKKS(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
      finalRes.correct = recon_pic.compare(pic);
      cout << "the invKT is: " << tools.invKT[0] << endl;
      if (writeSimpleData) {
        finalRes.writeResult(filename);
      }
    }
  }
  else if (mode == 4) {
    map<vector<int>, double> cal_ans;
    int bits;
    cout << "please input the repeat times and coeff_modulus bits:";
    cin >> bits;
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        vector<int> recovery_ids;
        vector<SharePic>UploadShares;
        for (int k = 0; k < 6; k++)
          if (k != i && k != j) {
            recovery_ids.push_back(k);
            UploadShares.push_back(shares[k]);
          }
        recoveryByPlainCKKS(UploadShares, picParms, tools);
        if (tools.invKT[0] > picParms.getP() / 2) {
          tools.invKT[0] -= picParms.getP();
        }
        picParms.setBaseLenOfModulus(bits);
        picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
        Picture recon_pic = recoryShareCKKS(picParms, UploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
        finalRes.correct = recon_pic.compare(pic);
        cout << "the invKT is: " << tools.invKT[0] << endl;
        cal_ans[recovery_ids] = finalRes.correct;
      }
    }
    fstream out("all.txt", ios::out);
    map<vector<int>, double>::iterator iter = cal_ans.begin();
    while (iter != cal_ans.end()) {
      vector<int> tem = iter->first;
      for (int i = 0; i < tem.size(); i++) {
        out << tem[i] << ",";
      }
      out << iter->second << endl;
      iter++;
    }
    out.close();
  }

  if (outPutCip)
  {
    Plaintext tem;
    vector<Ciphertext>x, fx;
    vector<Plaintext>x_plain, fx_plain;
    for (int i = 0; i < uploadShares.size(); i++) {
      x.push_back(uploadShares[i].x_en[0]);
      x_plain.push_back(uploadShares[i].X_plain);
      fx.push_back(uploadShares[i].fx_en[0]);
      fx_plain.push_back(uploadShares[i].fx_plain);
    }
    tem = x_plain[0];
    for (int i = 0; i < 10; i++) {
      cout << tem[i] << " ";
    }
    tem = fx_plain[0];

    for (int i = 0; i < 10; i++) {
      cout << tem[i] << " ";
    }
    cout << endl;
    cout << "output the shared pic:\n";
    outPutCiphertext(fx, "ciphertext");
    cout << "output the x_id:\n";
    outPutCiphertext(x, "xid");
    cout << "output the plainText of x\n";
    outPutCiphertext(x_plain, "x_plaintext");
    cout << "output the plainText of vector y";
    outPutCiphertext(fx_plain, "vector_fx");
  }
}

void CKKS2(int mode, Params& picParms, vector<SharePic>& uploadShares, DecTools& tools,
  Result& finalRes, Picture& pic, int writeSimpleData, string filename, ifstream& datain, vector<SharePic>& shares, int degree = 8192, int p_mo = 1024, bool printAns = true) {
  while (1) {
    cout << "please choose test mods:\n";
    cout << "------------------------------------\n";
    cout << "1.normal example\n" << "2.change the regular coeff modulus bits\n";
    cout << "3.repeat the normal example to test the ans's randomness:\n";
    cout << "4.check every possible ks' scheme:\n";
    cout << "0:exit\n";
    cout << "------------------------------------\n";
    cin >> mode;
    if (mode == 0)
      break;
    int outPutCip;
    cout << "whether to output the CiphereText(0|1):";
    cin >> outPutCip;
    if (mode == 1) {
      //每组x,y输入，然后就可以得到想要的y
  //recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
      picParms.setMaxLevel(2);
      picParms.setBaseLenOfModulus(40);
      picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
      Picture recon_pic = recoryShareCKKS2(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
      finalRes.correct = recon_pic.compare(pic);
      //cout << "the invKT is: " << tools.invKT[0] << endl;
      if (writeSimpleData) {
        finalRes.writeResult(filename);
      }
    }
    else if (mode == 2) {
      cout << "please input the range of coeff_modulus bits[start, step, end]";
      int start, step, end;
      cin >> start >> step >> end;
      for (int i = start; i <= end; i += step) {
        picParms.setBaseLenOfModulus(i);
        picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
        Picture recon_pic = recoryShareCKKS(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
        finalRes.correct = recon_pic.compare(pic);
        cout << "the invKT is: " << tools.invKT[0] << endl;
        if (writeSimpleData) {
          finalRes.writeResult(filename);
        }
      }
    }
    else if (mode == 3) {
      int ti, bits;
      cout << "please input the repeat times and coeff_modulus bits:";
      cin >> ti >> bits;
      while (ti--) {
        //每组x,y输入，然后就可以得到想要的y
        //recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
        picParms.setBaseLenOfModulus(bits);
        picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
        Picture recon_pic = recoryShareCKKS2(picParms, uploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
        finalRes.correct = recon_pic.compare(pic);
        cout << "the invKT is: " << tools.invKT[0] << endl;
        if (writeSimpleData) {
          finalRes.writeResult(filename);
        }
      }
    }
    else if (mode == 4) {
      map<vector<int>, double> cal_ans;
      int bits;
      cout << "please input the repeat times and coeff_modulus bits:";
      cin >> bits;
      for (int i = 0; i < 6; i++) {
        for (int j = i + 1; j < 6; j++) {
          vector<int> recovery_ids;
          vector<SharePic>UploadShares;
          for (int k = 0; k < 6; k++)
            if (k != i && k != j) {
              recovery_ids.push_back(k);
              UploadShares.push_back(shares[k]);
            }
          recoveryByPlainCKKS(UploadShares, picParms, tools);
          if (tools.invKT[0] > picParms.getP() / 2) {
            tools.invKT[0] -= picParms.getP();
          }
          picParms.setBaseLenOfModulus(bits);
          picParms.setScale(pow(2.0, picParms.getBaseLenOfModulus()));
          Picture recon_pic = recoryShareCKKS2(picParms, UploadShares, tools, finalRes, pic, finalRes.poly_d1, finalRes.plain_m1, printAns);
          finalRes.correct = recon_pic.compare(pic);
          cout << "the invKT is: " << tools.invKT[0] << endl;
          cal_ans[recovery_ids] = finalRes.correct;
        }
      }
      fstream out("all.txt", ios::out);
      map<vector<int>, double>::iterator iter = cal_ans.begin();
      while (iter != cal_ans.end()) {
        vector<int> tem = iter->first;
        for (int i = 0; i < tem.size(); i++) {
          out << tem[i] << ",";
        }
        out << iter->second << endl;
        iter++;
      }
      out.close();
    }
    if (outPutCip)
    {
      Plaintext tem;
      vector<Ciphertext>x, fx;
      vector<Plaintext>x_plain, fx_plain;
      for (int i = 0; i < uploadShares.size(); i++) {
        x.push_back(uploadShares[i].x_en[0]);
        x_plain.push_back(uploadShares[i].X_plain);
        fx.push_back(uploadShares[i].fx_en[0]);
        fx_plain.push_back(uploadShares[i].fx_plain);
      }
      tem = x_plain[0];
      for (int i = 0; i < 10; i++) {
        cout << tem[i] << " ";
      }
      tem = fx_plain[0];

      for (int i = 0; i < 10; i++) {
        cout << tem[i] << " ";
      }
      cout << endl;
      cout << "output the shared pic:\n";
      outPutCiphertext(fx, "ciphertext");
      cout << "output the x_id:\n";
      outPutCiphertext(x, "xid");
      cout << "output the plainText of x\n";
      outPutCiphertext(x_plain, "x_plaintext");
      cout << "output the plainText of vector y";
      outPutCiphertext(fx_plain, "vector_fx");
    }
  }
}


int main_Test() {
  vector<double> x, y, z,temAns;
  x = { 13478400.1496, -2.0, -3.0 };
  y = { 249, -3.0, -4.0 };
  z = { 1, -4.0, -5.0 };

  EncryptionParameters parms(scheme_type::CKKS);

  size_t poly_modulus_degree = 8192*4;
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 50, 50, 50, 50, 50, 50, 60 }));
  double scale = pow(2.0, 50);

  auto context = SEALContext::Create(parms);

  KeyGenerator keygen(context);
  auto public_key = keygen.public_key();
  auto secret_key = keygen.secret_key();
  auto relin_keys = keygen.relin_keys();

  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);

  CKKSEncoder encoder(context);
  Plaintext xp, yp, zp,temUse;
  encoder.encode(x, scale, xp);
  encoder.encode(y, scale, yp);
  encoder.encode(z, scale, zp);
  Ciphertext xc, yc, zc;
  encryptor.encrypt(xp, xc);
  encryptor.encrypt(yp, yc);
  encryptor.encrypt(zp, zc);
  Ciphertext temp;
  Ciphertext result_c;
  evaluator.multiply(xc, yc, temp);
  evaluator.relinearize_inplace(temp, relin_keys);
  evaluator.rescale_to_next_inplace(temp);
  decryptor.decrypt(temp, temUse);
  encoder.decode(temUse,temAns);
  for (int i = 0; i < min((int)temAns.size(), 3); i++)
    cout << temAns[i] << " ";
  cout << endl;

  Plaintext wt;
  encoder.encode(1.0, scale, wt);
  cout << "    + Modulus chain index for zc: "
    << context->get_context_data(zc.parms_id())->chain_index() << endl;
  cout << "    + Modulus chain index for temp(x*y): "
    << context->get_context_data(temp.parms_id())->chain_index() << endl;
  cout << "    + Modulus chain index for wt: "
    << context->get_context_data(wt.parms_id())->chain_index() << endl;

  evaluator.multiply_plain_inplace(zc, wt);
  evaluator.rescale_to_next_inplace(zc);

  decryptor.decrypt(zc, temUse);
  encoder.decode(temUse, temAns);
  cout << temAns.size() << endl;
  cout << encoder.slot_count() << endl;
  for (int i = 0; i < min((int)temAns.size(), 3); i++)
    cout << temAns[i] << "  ";
  cout << endl;

  cout << "    + Modulus chain index for zc after zc*wt and rescaling: "
    << context->get_context_data(zc.parms_id())->chain_index() << endl;

  evaluator.multiply_inplace(temp, zc);
  evaluator.relinearize_inplace(temp, relin_keys);
  evaluator.rescale_to_next_inplace(temp);


  Plaintext result_p;
  decryptor.decrypt(temp, result_p);
  vector<double> result;
  encoder.decode(result_p, result);
  cout << result.size()<<endl;
  for (int i = 0; i < min((int)result.size(),3); i++)
    cout << result[i] << "  ";
  cout << endl;
  return 0;
}

int main_test() {
  EncryptionParameters parms(scheme_type::CKKS);
  size_t poly_modulus_degree = 16384;
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 40, 40, 40, 40, 40, 60 }));

  double scale = pow(2.0, 40);

  shared_ptr<seal::SEALContext> context = SEALContext::Create(parms);

  KeyGenerator keygen(context);
  auto public_key = keygen.public_key();
  auto secret_key = keygen.secret_key();
  auto relin_keys = keygen.relin_keys();

  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);

  CKKSEncoder encoder(context);
  Plaintext TWO,ONE;
  Ciphertext t_en,one_en;
  encoder.encode(-1.3, scale, TWO);
  encryptor.encrypt(TWO, t_en);
  vector<Ciphertext> vec;
  for (int i = 0; i < 63; i++) {
    vec.push_back(t_en);
  }
  vector<Ciphertext> ONES;
  encoder.encode(1, scale, ONE);
  encryptor.encrypt(ONE, one_en);
  for (int i = 0; i <= 7; i++) {
    ONES.push_back(one_en);
  }
  t_en = cal_mut(vec, evaluator, relin_keys, context,ONES);
  decryptor.decrypt(t_en, TWO);
  vector<double> ans;
  encoder.decode(TWO, ans);
  cout << ans[1];

}

/// <summary>
/// the  main program of BFV scheme
/// </summary>
/// <returns></returns>
int main_bfv() {
  srand(time(0));

  ifstream datain("cin.txt", ios::in);
  //get the essential parameter
  int l,printAns,writeSimpleData;
  while (1) {
    cout << "please input the size(L,W) of pic and (t,n) and (P,x_range):";
    //cin >> l;
    datain >> l;
    if (l == 10000)
      break;
    Result finalRes;
    finalRes.l = l;
    cin >> finalRes.l >>finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    //datain >> finalRes.l >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    Params picParms(finalRes.l, finalRes.w, finalRes.p, finalRes.n, finalRes.t, finalRes.x_range);
    //get the base pic according to params
    cout << "please input the poly_modulus_degree and plain_modulus of getShare and reconstruction:";
    cin >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
    //datain >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;

    cout << "whether to print the detail pixes of pictures(0|1):";
    cin >> printAns;
    //datain >> printAns;
    cout << "whether to record the the simple ans of the ans (and filename): ";
    cin >> writeSimpleData;

    Picture pic(picParms);
    pic.generatePic(picParms);
    if (printAns) {
      cout << "the pixes of orign_pic are:";
      pic.printPic();
    }
    vector<SharePic>shares;
    int n = picParms.getN();
    for (int i = 0; i < n; i++) {
      SharePic picTem(picParms);
      shares.push_back(picTem);
    }
    getShareByPlain(pic, shares, picParms.getk(), picParms.getXrange(), printAns);

    /*suppose the max number of pixes are no bigger than 1000
    under the control of rangeX and randMul of orign Pic*/
    getShareByHE(finalRes.poly_d0, finalRes.plain_m0, pic, shares, picParms, finalRes, printAns);

    cout << endl << endl << endl;
    cout << "Please putin " << picParms.getk() << " indexs of shares to recovery the pictures" << endl;
    vector<int> recovery_ids;
    vector<SharePic>uploadShares;

    for (int i = 0; i < picParms.getk(); i++) {
      int tem;
      cin >> tem;
      //datain >> tem;
      recovery_ids.push_back(tem);
      uploadShares.push_back(shares[tem]);
    }

    DecTools tools;
    Picture pic1 = recoveryByPlain(uploadShares, picParms, tools);
    cout << endl << endl;
    cout << "the recovery result using plaintext are : ";
    if (printAns) {
      pic1.printPic();
    }
    pic1.compare(pic);

    cout << endl;
    //每组x,y输入，然后就可以得到想要的y
    //recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
    Picture recon_pic = recoryShare(picParms, uploadShares, tools, finalRes, finalRes.poly_d1, finalRes.plain_m1, printAns);
    finalRes.correct = recon_pic.compare(pic);
    if (writeSimpleData) {
      finalRes.writeResult();
    }
  }
  datain.close();
}

/// <summary>
/// the main program of CKKS scheme
/// </summary>
/// <returns></returns>
int main(){

  srand(time(0));
  string filename;
  ifstream datain("cin.txt", ios::in);
  //get the essential parameter
  int l,printAns,writeSimpleData;
  while (1) {
    cout << "please input the size(L,W) of pic and (t,n) and (P,x_range):";
    cin >> l;
    //datain >> l;
    if (l == 10000)
      break;
    Result finalRes;
    finalRes.l = l;
    cin >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    //datain >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    Params picParms(finalRes.l, finalRes.w, finalRes.p, finalRes.n, finalRes.t, finalRes.x_range);
    picParms.setMaxLevel(7); // the number of Coeff_Modulus;
    picParms.setBaseLenOfModulus(50);
    picParms.setScale(pow(2.0, 50));

    //get the base pic according to params
    cout << "please input the poly_modulus_degree and plain_modulus of getShare and reconstruction:";
    cin >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
    //datain >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
    cout << "whether to print the detail pixes of pictures(0|1):";
    cin >> printAns;
    //datain >> printAns;
    cout << "whether to record the the simple ans of the ans (and filename): ";
    cin >> writeSimpleData >> filename;

    Picture pic(picParms);
    pic.generatePic(picParms);
    if (printAns) {
      cout << "the pixes of orign_pic are:";
      pic.printPic();
    }
    vector<SharePic>shares;
    int n = picParms.getN();
    for (int i = 0; i < n; i++) {
      SharePic picTem(picParms);
      shares.push_back(picTem);
    }
    getShareByPlain_VEC(pic, shares, picParms.getk(), picParms.getXrange(), printAns);

 
    /*suppose the max number of pixes are no bigger than 1000
    under the control of rangeX and randMul of orign Pic*/
    getShareByCKKS(finalRes.poly_d0, pic, shares, picParms, finalRes, printAns);

    cout << endl << endl << endl;
    cout << "Please putin " << picParms.getk() << " indexs of shares to recovery the pictures" << endl;
    vector<int> recovery_ids;
    vector<SharePic>uploadShares;

    for (int i = 0; i < picParms.getk(); i++) {
      int tem;
      cin >> tem;
      //datain >> tem;
      recovery_ids.push_back(tem);
      uploadShares.push_back(shares[tem]);
    }

    //the X_ids are as follow:
   // cout << "the X_ids are as follow: ";
 ///  for (int i = 0; i < uploadShares.size(); i++) {
 //     cout << uploadShares[i].X[0] << " ";
 //   }
   // cout << endl;

    DecTools tools;
    Picture pic1 = recoveryByPlainCKKS(uploadShares, picParms, tools);
    cout << endl << endl;
    cout << "the recovery result by plaintext are : ";
    if (printAns) {
      pic1.printPic();
    }
    pic1.compare(pic);


    cout << endl;

    if (tools.invKT[0] > picParms.getP() / 2) {
      tools.invKT[0] -= picParms.getP();
    }
    while (1) {
      cout << "please choose test mods:\n";
      cout << "------------------------------------\n";
      cout << "1.normal example\n" << "2.change the regular coeff modulus bits\n";
      cout << "3.repeat the normal example to test the ans's randomness:\n";
      cout << "4.check every possible ks' scheme:\n";
      cout << "5.step into improvment state;\n";
      cout << "0:exit\n";
      cout << "------------------------------------\n";
      int tag;
      cin >> tag;
      if (tag == 0)
        break;
      if (tag == 5)
        CKKS2(tag, picParms, uploadShares, tools, finalRes, pic, writeSimpleData, filename, datain, shares, finalRes.poly_d1, finalRes.plain_m1, printAns);
      functionTest(tag, picParms, uploadShares, tools, finalRes, pic, writeSimpleData, filename, datain,shares, finalRes.poly_d1, finalRes.plain_m1, printAns);
    }
   }
  datain.close();
}

//
//int main_SIS2() {
//  srand(time(0));
//
//  //ifstream datain("cin.txt", ios::in);
//  //get the essential parameter
//  int l;
//  while (1) {
//    cout << "please input the size(L,W) of pic and (t,n) and (P,x_range):";
//    cin >> l;
//    //datain >> l;
//    if (l == 10000)
//      break;
//    Result finalRes;
//    finalRes.l = l;
//    cin >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
//    //datain >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
//    Params picParms(finalRes.l, finalRes.w, finalRes.p, finalRes.n, finalRes.t, finalRes.x_range);
//    picParms.setMaxLevel(7); // the number of Coeff_Modulus;
//    picParms.setBaseLenOfModulus(50);
//    picParms.setScale(pow(2.0, 50));
//
//    //get the base pic according to params
//    cout << "please input the poly_modulus_degree and plain_modulus of getShare and reconstruction:";
//    cin >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
//    //datain >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
//    Picture pic(picParms);
//    pic.generatePic(picParms);
//    cout << "the pixes of orign_pic are:";
//    pic.printPic();
//    vector<SharePic>shares;
//    int n = picParms.getN();
//    for (int i = 0; i < n; i++) {
//      SharePic picTem(picParms);
//      shares.push_back(picTem);
//    }
//    getShareByPlain_VEC(pic, shares, picParms.getk(), picParms.getXrange());
//
//
//    /*suppose the max number of pixes are no bigger than 1000
//    under the control of rangeX and randMul of orign Pic*/
//    getShareByCKKS(finalRes.poly_d0, pic, shares, picParms, finalRes);
//
//    cout << endl << endl << endl;
//    cout << "Please putin " << picParms.getk() << " indexs of shares to recovery the pictures" << endl;
//    vector<int> recovery_ids;
//    vector<SharePic>uploadShares;
//
//    for (int i = 0; i < picParms.getk(); i++) {
//      int tem;
//      cin >> tem;
//      //datain >> tem;
//      recovery_ids.push_back(tem);
//      uploadShares.push_back(shares[tem]);
//    }
//
//    //the X_ids are as follow:
//   // cout << "the X_ids are as follow: ";
//    for (int i = 0; i < uploadShares.size(); i++) {
//      cout << uploadShares[i].X[0] << " ";
//    }
//    cout << endl;
//
//    DecTools tools;
//    Picture pic1 = recoveryByPlainCKKS(uploadShares, picParms, tools);
//    cout << endl << endl;
//    cout << "the recovery result are : ";
//    pic1.printPic();
//    pic1.compare(pic);
//
//
//    cout << endl;
//    //每组x,y输入，然后就可以得到想要的y
//    //recovery(X, Ys, m, K, KT, invKT, 8192, 1024);
//    Picture recon_pic = recoryShareCKKS(picParms, uploadShares, tools, finalRes, finalRes.poly_d1, finalRes.plain_m1);
//    finalRes.correct = recon_pic.compare(pic);
//    //finalRes.writeResult("CKKS.text");
//  }
//  ////datain.close();
//}
int main_BFV() {
  srand(time(0));

  string filename;
  ifstream datain("cin.txt", ios::in);
  //get the essential parameter
  int l, printAns, writeSimpleData;
  while (1) {
    cout << "please input the size(L,W) of pic and (t,n) and (P,x_range):";
    cin >> l;
    //datain >> l;
    if (l == 10000)
      break;
    Result finalRes;
    finalRes.l = l;
    cin >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    //datain >> finalRes.w >> finalRes.t >> finalRes.n >> finalRes.p >> finalRes.x_range;
    Params picParms(finalRes.l, finalRes.w, finalRes.p, finalRes.n, finalRes.t, finalRes.x_range);
    picParms.setMaxLevel(7); // the number of Coeff_Modulus;
    picParms.setBaseLenOfModulus(50);
    picParms.setScale(pow(2.0, 50));

    //get the base pic according to params
    cout << "please input the poly_modulus_degree and plain_modulus of getShare and reconstruction:";
    cin >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
    //datain >> finalRes.poly_d0 >> finalRes.plain_m0 >> finalRes.poly_d1 >> finalRes.plain_m1;
    cout << "whether to print the detail pixes of pictures(0|1):";
    cin >> printAns;
    //datain >> printAns;
    cout << "whether to record the the simple ans of the ans (and filename): ";
    cin >> writeSimpleData >> filename;

    Picture pic(picParms);
    pic.generatePic(picParms);
    if (printAns) {
      cout << "the pixes of orign_pic are:";
      pic.printPic();
    }
    vector<SharePic>shares;
    int n = picParms.getN();
    for (int i = 0; i < n; i++) {
      SharePic picTem(picParms);
      shares.push_back(picTem);
    }
    getShareByPlain_VEC(pic, shares, picParms.getk(), picParms.getXrange(), printAns);


    /*suppose the max number of pixes are no bigger than 1000
    under the control of rangeX and randMul of orign Pic*/
    getShareByBFV(finalRes.poly_d0, pic, shares, picParms, finalRes, printAns);

    cout << endl << endl << endl;
    cout << "Please putin " << picParms.getk() << " indexs of shares to recovery the pictures" << endl;
    vector<int> recovery_ids;
    vector<SharePic>uploadShares;

    for (int i = 0; i < picParms.getk(); i++) {
      int tem;
      cin >> tem;
      //datain >> tem;
      recovery_ids.push_back(tem);
      uploadShares.push_back(shares[tem]);
    }

    //the X_ids are as follow:
   // cout << "the X_ids are as follow: ";
 ///  for (int i = 0; i < uploadShares.size(); i++) {
 //     cout << uploadShares[i].X[0] << " ";
 //   }
   // cout << endl;

    DecTools tools;
    Picture pic1 = recoveryByPlainCKKS(uploadShares, picParms, tools);
    cout << endl << endl;
    cout << "the recovery result by plaintext are : ";
    if (printAns) {
      pic1.printPic();
    }
    pic1.compare(pic);


    cout << endl;

    if (tools.invKT[0] > picParms.getP() / 2) {
      tools.invKT[0] -= picParms.getP();
    }
    while (1) {
      cout << "please choose test mods:\n";
      cout << "------------------------------------\n";
      cout << "1.normal example\n" << "2.change the regular coeff modulus bits\n";
      cout << "3.repeat the normal example to test the ans's randomness:\n";
      cout << "4.check every possible ks' scheme:\n";
      cout << "0:exit\n";
      cout << "------------------------------------\n";
      int tag;
      cin >> tag;
      if (tag == 0)
        break;
      functionTestBFV(tag, picParms, uploadShares, tools, finalRes, pic, writeSimpleData, filename, datain, shares, finalRes.poly_d1, finalRes.plain_m1, printAns);
    }
  }
  datain.close();
}
