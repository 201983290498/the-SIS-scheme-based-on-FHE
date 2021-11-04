#pragma once
#include"baseClass.h"
using namespace std;

/**
  the default poly_degree is 4
*/
void getShareByPlain(Picture& source, vector<SharePic>& shares, int k = 4, int range_x = 10, bool printAns = true) {
  vector<ll> y = source.getSec();
  int n = shares.size();

  // use poly format to generate the sharePic
  for (int index = 0; index < y.size(); index = index + k) {
    vector<ll> temy, temx;
    for (int i = 0; i < k; i++)
      temy.push_back(y[index + i]);
    //generate the shared Pixes in sharePic
    //escape existing tow same number
    for (int pn = 0; pn < n; pn++) {
      int x = rand() % range_x;
      while (1) {
        bool flag = true;
        for (int i = 0; i < temx.size(); i++) {
          if (temx[i] == x)
            flag = false;
        }
        if (flag)
          break;
        else
          x = rand() % range_x;
      }
      shares[pn].addNewPixByPlain(temy, x);
      temx.push_back(x);
    }
  }
  if (printAns) {
    for (int i = 0; i < n; i++) {
      cout << "the " << i << "shares are:";
      shares[i].PrintFx();
    }
  }
}

// they have common X
void getShareByPlain_VEC(Picture& source, vector<SharePic>& shares, int k = 4, int range_x = 10, bool printAns = true) {
  vector<ll> y = source.getSec();
  vector<int> X;
  int n = shares.size();

  //set the common X of the shares and
  // the hidden information length of X; 
  for (int pn = 0; pn < n; pn++) {
    int x = rand() % range_x;
    while (1) {
      bool flag = true;
      for (int i = 0; i < X.size(); i++) {
        if (X[i] == x)
          flag = false;
      }
      if (flag)
        break;
      else
        x = rand() % range_x;
    }
    X.push_back(x);
    shares[pn].X.push_back(x);
    // the size of X_len is to overcome the conflication between y.size and k
    shares[pn].X_len = (y.size() % k == 0 ? y.size() / k : (y.size() - y.size() % k + k) / k);
  }
 

  // use poly format to generate the sharePic
  for (int index = 0; index < y.size(); index = index + k) {
    vector<ll> temy;
    for (int i = 0; i < k; i++)
      temy.push_back(y[index + i]);
    //generate the shared Pixes in sharePic
    //escape existing tow same number
    for (int pn = 0; pn < n; pn++)
      shares[pn].addNewPixByPlain(temy,X[pn],false);
  }
  if (printAns) {
    for (int i = 0; i < n; i++) {
      cout << "the " << i << "shares are:";
      shares[i].PrintFx();
    }
  }
}


vector<int> recoveryDemo(vector<int>& Xtem, vector<int>& ytem, DecTools& tools, Params& parms, int index) {
  int tem;
  ll P = parms.getP(), k = parms.getk();
  if (index < 64 * 64) {
    ll& KT = tools.KT[index];
    vector<ll>& m = tools.m[index]; m.clear();
    vector<ll>& K = tools.K[index]; K.clear();
    vector<ll>& invm = tools.invm[index]; invm.clear();
    vector<ll>& A = tools.A[index]; A.clear();
    ll& invKT = tools.invKT[index]; 
    KT = 1;
    for (int i = 0; i < parms.getk(); i++) {//cacul m ，k
      tem = 1;
      int xi = Xtem[i];
      for (int j = 0; j < parms.getk(); j++) {
        if (j == i) {
          continue;
        }
        int xj = Xtem[j];
        tem *= (xi - xj);
      }
      m.push_back(tem);
      /// <summary>
      /// print M
      /// </summary>
      /// <param name="Xtem"></param>
      /// <param name="ytem"></param>
      /// <param name="tools"></param>
      /// <param name="parms"></param>
      /// <param name="index"></param>
      /// <returns></returns>

     // cout << "m:" << tem<<" ";
      KT = KT * tem % parms.getP(); //cacul K
    }
    // cout << endl;
     /// <summary>
     /// print A
     /// </summary>
     /// <param name="Xtem"></param>
     /// <param name="ytem"></param>
     /// <param name="tools"></param>
     /// <param name="parms"></param>
     /// <param name="index"></param>
     /// <returns></returns>
    for (int i = 0; i < parms.getk(); i++) {
      long long temm = mypow(m[i], P - 2, P);//m[i] can be a nagetive
      long long temK = KT * temm % P;
      K.push_back(temK);//push in Ki
      invm.push_back(temm);//push in imK,it can be a negative
      A.push_back(temK * ytem[i]);//cal
      //cout << (temK * ytem[i]) << " ";
    }
    // cout << endl;
    vector<int> D = getD(Xtem, k);//here is a problem
    vector<int> a;
    invKT = mypow(KT, P - 2, P);//the result is mod P
    /// <summary>
    /// print invKT
    /// </summary>
    /// <param name="Xtem"></param>
    /// <param name="ytem"></param>
    /// <param name="tools"></param>
    /// <param name="parms"></param>
    /// <param name="index"></param>
    /// <returns></returns>
   // cout << "invKT is " << invKT << endl;
    for (int i = 0; i < k; i++) {
      int tem = 0;
      for (int j = 0; j < k; j++) {
        tem += A[j] * D[j * k + k - 1 - i] % P;//the result maybe vary large
      }
      //cout << tem << " ";
      tem = tem * invKT;
      a.push_back((tem % P + P) % P);
    }
    //cout << endl;
    return a;
  }
  else {
    ll KT;
    vector<ll> m;
    vector<ll> K;
    vector<ll> invm;
    vector<ll> A;
    ll invKT;
    KT = 1;
    for (int i = 0; i < parms.getk(); i++) {//cacul m ，k
      tem = 1;
      int xi = Xtem[i];
      for (int j = 0; j < parms.getk(); j++) {
        if (j == i) {
          continue;
        }
        int xj = Xtem[j];
        tem *= (xi - xj);
      }
      m.push_back(tem);
      /// <summary>
      /// print M
      /// </summary>
      /// <param name="Xtem"></param>
      /// <param name="ytem"></param>
      /// <param name="tools"></param>
      /// <param name="parms"></param>
      /// <param name="index"></param>
      /// <returns></returns>

     // cout << "m:" << tem<<" ";
      KT = KT * tem % parms.getP(); //cacul K
    }
    // cout << endl;
     /// <summary>
     /// print A
     /// </summary>
     /// <param name="Xtem"></param>
     /// <param name="ytem"></param>
     /// <param name="tools"></param>
     /// <param name="parms"></param>
     /// <param name="index"></param>
     /// <returns></returns>
    for (int i = 0; i < parms.getk(); i++) {
      long long temm = mypow(m[i], P - 2, P);//m[i] can be a nagetive
      long long temK = KT * temm % P;
      K.push_back(temK);//push in Ki
      invm.push_back(temm);//push in imK,it can be a negative
      A.push_back(temK * ytem[i]);//cal
      //cout << (temK * ytem[i]) << " ";
    }
    // cout << endl;
    vector<int> D = getD(Xtem, k);//here is a problem
    vector<int> a;
    invKT = mypow(KT, P - 2, P);//the result is mod P
    /// <summary>
    /// print invKT
    /// </summary>
    /// <param name="Xtem"></param>
    /// <param name="ytem"></param>
    /// <param name="tools"></param>
    /// <param name="parms"></param>
    /// <param name="index"></param>
    /// <returns></returns>
   // cout << "invKT is " << invKT << endl;
    for (int i = 0; i < k; i++) {
      int tem = 0;
      for (int j = 0; j < k; j++) {
        tem += A[j] * D[j * k + k - 1 - i] % P;//the result maybe vary large
      }
      //cout << tem << " ";
      tem = tem * invKT;
      a.push_back((tem % P + P) % P);
    }
    //cout << endl;
    return a;
  }
}

Picture& recoveryByPlain(vector<SharePic>& source, Params& parms, DecTools& tools) {
  cout << "----------------";
  cout << "----------------";
  cout << "start recover the secret image using plaintext without FHE:";
  int pixIndex = 0;
  Picture pic(parms);
  vector<int>Xtem, Ytem;
  for (int i = 0; i < source[0].X.size(); i++) {
    Xtem.clear(); Ytem.clear();
    for (int j = 0; j < source.size(); j++) {
      Xtem.push_back(source[j].X[i]);
      Ytem.push_back(source[j].fx[i]);
    }
    vector<int> pixes = recoveryDemo(Xtem, Ytem, tools, parms, i);
    for (int j = 0; j < pixes.size(); j++) {
      pic.pushPies(pixes[j], pixIndex++);
    }
  }
  return pic;
}


Picture recoveryByPlainCKKS(vector<SharePic>& source, Params& parms, DecTools& tools) {
  cout << "----------------\n";
  cout << "----------------\n";
  cout << "start recover the secret image using plaintext without FHE:";
  int pixIndex = 0;
  Picture pic(parms);
  vector<int>Xtem, Ytem;
  for (int i = 0; i < source[0].fx.size(); i++) {
    Xtem.clear(); Ytem.clear();
    for (int j = 0; j < source.size(); j++) {
      Xtem.push_back(source[j].X[0]);
      Ytem.push_back(source[j].fx[i]);
    }
    vector<int> pixes = recoveryDemo(Xtem, Ytem, tools, parms, i);
    for (int j = 0; j < pixes.size(); j++) {
      pic.pushPies(pixes[j], pixIndex++);
    }
  }
  return pic;
}

