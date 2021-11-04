#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include"tools.h"
#include "seal/seal.h"
using namespace std;
using namespace seal;

using parms_id_type = util::HashFunction::sha3_block_type;

// the use of context is to check the chain-index of Ciphertext;
int limit[] = { 2, 4, 8, 16, 32, 64, 128 };

/// <summary>
/// get the level of Ciplhertext
/// </summary>
/// <param name="tem"></param>
/// <param name="context"></param>
/// <param name="tag">
/// if tag=0 ,it means add, which is default situation;when tag = 1, it means mul
/// </param>
/// <returns></returns>
int getLevel(Ciphertext& tem, shared_ptr<seal::SEALContext>& context, int tag=0) {
  int level = context->get_context_data(tem.parms_id())->chain_index();
  if (tag && !level) {
    cout << "the scale is not support multiply again";
    throw "the scale is not support multiply again";
  }
  else {
    return level;
  }
}

void modifyCipScale(Evaluator& evaluator, Ciphertext& tem,parms_id_type param_id,double scale) {
  evaluator.mod_switch_to_inplace(tem, param_id);
  tem.scale() = scale;
}

int reSize(int tem) {
  int i;
  for (i = 0; i < 7; i++) {
    if (tem <= limit[i])
      break;
  }
  return limit[i];
}

void multiply(Ciphertext& tem1, Ciphertext tem2, Evaluator& evaluator, RelinKeys& relin_key) {
  ios old_fmt(nullptr);
  old_fmt.copyfmt(cout);
  cout << fixed << setprecision(10);
  evaluator.multiply_inplace(tem1, tem2);
  evaluator.relinearize_inplace(tem1, relin_key);
  evaluator.rescale_to_next_inplace(tem1);
  cout.copyfmt(old_fmt);
}

void multiply1(Ciphertext& tem1, Ciphertext tem2, Evaluator& evaluator, RelinKeys& relin_key) {
  evaluator.multiply_inplace(tem1, tem2);
  evaluator.relinearize_inplace(tem1, relin_key);
}


//mantually change the scale of Ciphertext and using mod_switch to change level 
void fitSize2(vector<Ciphertext>& vec, shared_ptr<seal::SEALContext>& context, vector<Ciphertext>& ONES, Evaluator& evaluator, RelinKeys& relin_keys) {
  if (vec.size() == 1)
    return;
  // change the scale and chain-index of  every element
  int minn = 10, tem, minn_id = 0;
  parms_id_type minn_p;
  for (int i = 0; i < vec.size(); i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    //change point 1
    if (minn > tem) {
      minn = tem;
      minn_p = vec[i].parms_id();
      minn_id = i;
    }
  }
  int  siz = vec.size();
  int real_siz = reSize(siz);

  // if the chain_index is not the lowest index,
  //we should change the chain-index of Ciphertext
  for (int i = 0; i < siz; i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    if (tem > minn) {
      // or using 2^40;
      modifyCipScale(evaluator, vec[i], minn_p, vec[minn_id].scale());
    }
  }
  for (int i = siz; i < real_siz; i++) {
    vec.push_back(ONES[minn]);
  }
}


// the ONES is used to change Ciphertext.scale and reSize the vector
void fitSize(vector<Ciphertext>& vec, shared_ptr<seal::SEALContext>& context, vector<Ciphertext>& ONES, Evaluator& evaluator, RelinKeys& relin_keys) {
  if (vec.size() == 1)
    return;
  // change the scale and chain-index of  every element
  int minn = 10, tem;
  parms_id_type minn_p;
  for (int i = 0; i < vec.size(); i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    minn = minn < tem ? minn : tem;
  }
  int  siz = vec.size();
  int real_siz = reSize(siz);

  // if the chain_index is not the lowest index,we should change the chain-index of Ciphertext
  for (int i = 0; i < siz; i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    while (tem > minn) {
      multiply(vec[i], ONES[tem], evaluator, relin_keys);
      tem--;
    }
  }
  for (int i = siz; i < real_siz; i++) {
    vec.push_back(ONES[minn]);
  }
}

Ciphertext  cal_mut(vector<Ciphertext> vec, Evaluator& evaluator, RelinKeys& relin_keys, shared_ptr<seal::SEALContext>& context,vector<Ciphertext>& ONES) {
  vector<Ciphertext> pre_ans;
  // here is a problem: the ONES has to be contructed earlier; 
  fitSize(vec, context, ONES, evaluator, relin_keys);
  while (vec.size() > 1) {
    for (int i = 0; i < vec.size(); i += 2) {
      Ciphertext tem = vec[i];
      //cout << context->get_context_data(tem.parms_id())->chain_index() << " " << tem.scale() << "    ";
      multiply(tem, vec[i + 1], evaluator, relin_keys);
      pre_ans.push_back(tem);
      //cout << context->get_context_data(tem.parms_id())->chain_index() << " " << tem.scale() << endl;
    }
    vec = pre_ans;
    pre_ans.clear();
  }
  return vec[0];
}

/// <summary>
/// the add funciotn of CKKS
/// </summary>
/// <param name="evaluator"></param>
/// <param name="destinction"></param>
/// <param name="tem2"></param>
/// <param name="context"></param>
void add(Evaluator& evaluator, Ciphertext& destinction, Ciphertext& tem2, shared_ptr<seal::SEALContext>& context) {
  int level1 = getLevel(destinction, context);
  int level2 = getLevel(tem2, context);
  if (level1 < level2) {
    modifyCipScale(evaluator, tem2, destinction.parms_id(), destinction.scale());
  }
  else if(level1>level2) {
    modifyCipScale(evaluator, destinction, tem2.parms_id(), tem2.scale());
  }
  evaluator.add_inplace(destinction, tem2);
}

void sub(Evaluator& evaluator, Ciphertext& destinction, Ciphertext& tem2, shared_ptr<seal::SEALContext>& context) {
  int level1 = getLevel(destinction, context);
  int level2 = getLevel(tem2, context);
  if (level1 < level2) {
    modifyCipScale(evaluator, tem2, destinction.parms_id(), destinction.scale());
  }
  else if (level1 > level2) {
    modifyCipScale(evaluator, destinction, tem2.parms_id(), tem2.scale());
  }
  evaluator.sub_inplace(destinction, tem2);
}

void add(Evaluator& evaluator, Ciphertext& tem1, Ciphertext& tem2, Ciphertext& destinction, shared_ptr<seal::SEALContext>& context) {
  int level1 = getLevel(tem1, context);
  int level2 = getLevel(tem2, context);
  if (level1 < level2) {
    modifyCipScale(evaluator, tem2, tem1.parms_id(), tem1.scale());
  }
  else if (level1 > level2) {
    modifyCipScale(evaluator, tem1, tem2.parms_id(), tem2.scale());
  }
  evaluator.add(tem1, tem2,destinction);
}

void sub(Evaluator& evaluator, Ciphertext& tem1, Ciphertext& tem2, Ciphertext& destinction, shared_ptr<seal::SEALContext>& context) {
  int level1 = getLevel(tem1, context);
  int level2 = getLevel(tem2, context);
  if (level1 < level2) {
    modifyCipScale(evaluator, tem2, tem1.parms_id(), tem1.scale());
  }
  else if (level1 > level2) {
    modifyCipScale(evaluator, tem1, tem2.parms_id(), tem2.scale());
  }
  evaluator.sub(tem1, tem2, destinction);
}

int around(double x) {
  return (x > 0 ? floor(x + 0.5) : ceil(x - 0.5));
}

// the ONES is used to change Ciphertext.scale and reSize the vector
void fitSize(vector<Ciphertext>& vec, Ciphertext& ONES) {
  if (vec.size() == 1)
    return;

  // modify the number of vec to be 2^n
  int  siz = vec.size();
  int real_siz = reSize(siz);
  for (int i = siz; i < real_siz; i++) {
    vec.push_back(ONES);
  }
}

Ciphertext  cal_mutBFV(vector<Ciphertext> vec, Evaluator& evaluator, RelinKeys& relin_keys, Ciphertext& ONES) {
  vector<Ciphertext> pre_ans;
  // here is a problem: the ONES has to be contructed earlier; 
  fitSize(vec, ONES);
  while (vec.size() > 1) {
    for (int i = 0; i < vec.size(); i += 2) {
      Ciphertext tem = vec[i];
      multiply1(tem, vec[i + 1], evaluator, relin_keys);
      pre_ans.push_back(tem);
    }
    vec = pre_ans;
    pre_ans.clear();
  }
  return vec[0];
}


// the ONES is used to change Ciphertext.scale and reSize the vector
void fitSize4(vector<Ciphertext>& vec, shared_ptr<seal::SEALContext>& context, vector<Ciphertext>& ONES, Evaluator& evaluator, RelinKeys& relin_keys) {
  if (vec.size() == 1)
    return;
  // change the scale and chain-index of  every element
  int minn = 10, tem;
  parms_id_type minn_p;
  for (int i = 0; i < vec.size(); i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    minn = minn < tem ? minn : tem;
  }
  int  siz = vec.size();
  int real_siz = reSize(siz);

  // if the chain_index is not the lowest index,we should change the chain-index of Ciphertext
  for (int i = 0; i < siz; i++) {
    tem = context->get_context_data(vec[i].parms_id())->chain_index();
    while (tem > minn) {
      //cout << context->get_context_data(vec[i].parms_id())->chain_index() << " " << context->get_context_data(ONES[tem].parms_id())->chain_index();
      multiply(vec[i], ONES[tem], evaluator, relin_keys);
      tem--;
    }
  }
  for (int i = siz; i < real_siz; i++) {
    vec.push_back(ONES[minn]);
  }
}

Ciphertext  cal_mut1(vector<Ciphertext> vec, Evaluator& evaluator, RelinKeys& relin_keys, shared_ptr<seal::SEALContext>& context, vector<Ciphertext>& ONES, Decryptor& decryptor, CKKSEncoder& encoder) {
  vector<Ciphertext> pre_ans;
  // here is a problem: the ONES has to be contructed earlier; 
  fitSize4(vec, context, ONES, evaluator, relin_keys);
  while (vec.size() > 1) {
    for (int i = 0; i < vec.size(); i += 2) {
      Ciphertext tem = vec[i];
      //cout << context->get_context_data(tem.parms_id())->chain_index() << " " << tem.scale() << "    ";
      multiply(tem, vec[i + 1], evaluator, relin_keys);
      pre_ans.push_back(tem);
      //PrintDeVec(decryptor, encoder, tem);
      //cout << context->get_context_data(tem.parms_id())->chain_index() << " " << tem.scale() << endl;
    }
    vec = pre_ans;
    pre_ans.clear();
  }
  return vec[0];
}


