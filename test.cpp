#include <iostream>
#include <cassert>
#include <time.h>

#include "params.h"
#include "lwehe.h"

#include <iostream>
#include <fstream>
using namespace std;


void test_nand(SchemeLWE& s){
    Ctxt_LWE ct_res, ct_res01, ct_res23, ct0, ct1;
    KeyGen keygen(parLWE);
    SKey_base_LWE sk_base; //s
    SKey_boot sk_boot; //f
    BSKey_LWE bsk;
    HybridKey hybridkey, hybridkey_b;
    SKey_ks sk_ks; //z
    KSKey_LWE_vec ksk;
    NGSFFTctxt crs_a, crs_a_neg;

    int parties = 2;
    int for_num = 50;
    cout << "---- " << parties << " parties ----" << endl;

    vector<KSKey_LWE_vec> ksk_all;
    vector<BSKey_LWE> bsk_all;
    vector<HybridKey> hybridkey_all;
    vector<SKey_base_LWE> sk_base_all;

    cout << "Begin KeyGen--" << endl;

    keygen.get_crs_a(crs_a, crs_a_neg);
    hybridkey_b.b = crs_a_neg;

    float avg_time_kg = 0.0;
    for (size_t i = 0; i < parties; i++)
    {
        clock_t start_kg = clock();
        keygen.get_sk_base(sk_base); //s
        keygen.get_sk_ks(sk_ks); //z
        keygen.get_ksk_mk(ksk, sk_base, sk_ks);
        keygen.get_sk_boot(sk_boot);
        keygen.get_bsk_mk_binary(bsk, sk_base, sk_boot);
        keygen.get_hybridkey(hybridkey, sk_ks, sk_boot, crs_a);
        avg_time_kg += float(clock()-start_kg)/CLOCKS_PER_SEC;
        ksk_all.push_back(ksk);
        bsk_all.push_back(bsk);
        hybridkey_all.push_back(hybridkey);
        sk_base_all.push_back(sk_base);
    }

    hybridkey_all.push_back(hybridkey_b);

    cout << "KG time (per party): " << avg_time_kg/parties << endl;

    cout << "Begin NAND Test--" << endl;

    float avg_time = 0.0;
    int err_num = 0;
    for (int i = 0; i < for_num; i++){

        int m0 = binary_sampler(rand_engine);
        int m1 = binary_sampler(rand_engine);
        s.encrypt_mk(ct0, m0, sk_base_all);
        s.encrypt_mk(ct1, m1, sk_base_all);
        clock_t start01 = clock();
        s.nand_gate(ct_res01, ct0, ct1);
        s.mk_bootstrap_binary(ct_res01, bsk_all, hybridkey_all, ksk_all);
        avg_time += float(clock()-start01)/CLOCKS_PER_SEC;

        int m01 = s.decrypt(ct_res01, sk_base_all, Param::t);
        if (m01 != !(m0 & m1))
        {
            err_num ++;
        }
    }
    cout << "NAND+BS time: " << avg_time/for_num << endl;
    cout << "Error rate: " << float(err_num)/for_num << endl;

}


int main(){
    Param param(LWE);
    SchemeLWE s;
    test_nand(s);
}