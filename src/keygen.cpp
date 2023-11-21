#include "keygen.h"
#include "params.h"
#include "fft.h"

#include<iostream>
#include<time.h>
#include<algorithm>
#include<math.h>

using namespace NTL;
using namespace std;

void KeyGen::get_crs_a(NGSFFTctxt& crs_a, NGSFFTctxt& crs_a_neg){
    crs_a = NGSFFTctxt(param.l_hybrid, FFTPoly(param.N2p1));
    crs_a_neg = NGSFFTctxt(param.l_hybrid, FFTPoly(param.N2p1));
    vector<int> crs_a_tmp(param.N,0L);
    vector<int> crs_a_neg_tmp(param.N,0L);
    for (size_t i = 0; i < param.l_hybrid; i++)
    {
        sampler.get_uniform_vector_with_neg(crs_a_tmp, crs_a_neg_tmp);
        fftN.to_fft(crs_a[i], crs_a_tmp);
        fftN.to_fft(crs_a_neg[i], crs_a_neg_tmp);
    }

}

void KeyGen::get_sk_hybrid(SKey_hybrid& sk_hybrid){
    sk_hybrid = ModQPoly(Param::N,0);
    sampler.get_ternary_hwt_vector(sk_hybrid, param.h);
}

void KeyGen::get_hybridkey(HybridKey& hybridkey, SKey_hybrid& sk_hybrid, SKey_boot& sk_boot, NGSFFTctxt& crs_a){

    normal_distribution<double> gaussian_sampler(0.0, param.hybrid_e_st_dev);
    int tmp_e;
    FFTPoly sk_hybrid_fft(Param::N2p1);
    fftN.to_fft(sk_hybrid_fft, sk_hybrid);

    // generate hybridkey.b
    hybridkey.b = NGSFFTctxt(param.l_hybrid, FFTPoly(Param::N2p1));
    FFTPoly tmp_b_fft(Param::N2p1);
    vector<long> tmp_b(Param::N);
    vector<int> tmp_b_int(Param::N);
    for (size_t i = 0; i < param.l_hybrid; i++)
    {
        tmp_b_fft = crs_a[i] * sk_hybrid_fft;
        fftN.from_fft(tmp_b, tmp_b_fft);
        mod_q_boot(tmp_b_int, tmp_b);
        for (size_t j = 0; j < Param::N; j++)
        {
            tmp_e = static_cast<int>(round(gaussian_sampler(rand_engine)));
            tmp_b_int[j] = tmp_e - tmp_b_int[j];
        }
        fftN.to_fft(hybridkey.b[i], tmp_b_int);
    }

    // generate hybridkey.d_1
    hybridkey.d_1 = NGSFFTctxt(param.l_hybrid, FFTPoly(Param::N2p1));
    vector<int> tmp_d_1(Param::N);
    for (size_t i = 0; i < param.l_hybrid; i++){
        sampler.get_uniform_vector(tmp_d_1);
        fftN.to_fft(hybridkey.d_1[i], tmp_d_1);
    }

    // generate r
    SKey_hybrid r_hybrid(Param::N,0);
    sampler.get_ternary_hwt_vector(r_hybrid, param.h_hybrid);

    // generate hybridkey.d_0
    hybridkey.d_0 = NGSFFTctxt(param.l_hybrid, FFTPoly(Param::N2p1));
    FFTPoly tmp_d_0_fft(Param::N2p1);
    vector<long> tmp_d_0(Param::N);
    vector<int> tmp_d_0_int(Param::N);
    for (size_t i = 0; i < param.l_hybrid; i++)
    {
        tmp_d_0_fft = hybridkey.d_1[i] * sk_hybrid_fft;
        fftN.from_fft(tmp_d_0, tmp_d_0_fft);
        mod_q_boot(tmp_d_0_int, tmp_d_0);
        for (size_t j = 0; j < Param::N; j++)
        {
            tmp_e = static_cast<int>(round(gaussian_sampler(rand_engine)));
            tmp_d_0[j] = tmp_e - tmp_d_0_int[j] + pow(param.B_hybrid, i) * r_hybrid[j];
        }
        mod_q_boot(tmp_d_0_int, tmp_d_0);
        fftN.to_fft(hybridkey.d_0[i], tmp_d_0_int);
    }

    // generate hybridkey.d_2
    hybridkey.d_2 = NGSFFTctxt(param.l_hybrid, FFTPoly(Param::N2p1));
    FFTPoly tmp_d_2_fft(Param::N2p1);
    FFTPoly r_fft(Param::N2p1);
    vector<long> tmp_d_2(Param::N);
    vector<int> tmp_d_2_int(Param::N);
    for (size_t i = 0; i < param.l_hybrid; i++)
    {
        fftN.to_fft(r_fft, r_hybrid);
        tmp_d_2_fft = r_fft * crs_a[i];
        fftN.from_fft(tmp_d_2, tmp_d_2_fft);
        mod_q_boot(tmp_d_2_int, tmp_d_2);
        for (size_t j = 0; j < Param::N; j++)
        {
            tmp_e = static_cast<int>(round(gaussian_sampler(rand_engine)));
            tmp_d_2[j] = tmp_e + tmp_d_2_int[j] + pow(param.B_hybrid, i) * sk_boot.sk[j];
        }
        mod_q_boot(tmp_d_2_int, tmp_d_2);
        fftN.to_fft(hybridkey.d_2[i], tmp_d_2_int);
    }
}
void KeyGen::get_sk_boot(SKey_boot& sk_boot)
{
    // cout << "Started generating the secret key of the bootstrapping scheme" << endl;
    // clock_t start = clock();
    sk_boot.sk = ModQPoly(Param::N,0);
    sk_boot.sk_inv = ModQPoly(Param::N,0);

    sampler.get_invertible_vector(sk_boot.sk, sk_boot.sk_inv, Param::t, 1L);

    // cout << "sk_boot.sk: ";
    // for (size_t i = 0; i < sk_boot.sk.size(); i++)
    // {
    //     cout << sk_boot.sk[i] << " ";
    // }
    // cout << ". " << endl;

    // cout << "sk_boot.sk_inv: ";
    // for (size_t i = 0; i < sk_boot.sk_inv.size(); i++)
    // {
    //     cout << sk_boot.sk_inv[i] << " ";
    // }
    // cout << ". " << endl;

    // cout << "Generation time of the secret key of the bootstrapping scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_sk_ks(SKey_ks& sk_ks)
{
    // cout << "Started generating the secret key of the base scheme" << endl;
    // clock_t start = clock();
    sk_ks.clear();

    sk_ks = vector<int>(param.N, 0L);
    sampler.get_ternary_hwt_vector(sk_ks, param.h_hybrid);
    // cout << "Generation time of the secret key of the base scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_sk_base(SKey_base_LWE& sk_base)
{
    // cout << "Started generating the secret key of the base scheme" << endl;
    // clock_t start = clock();
    sk_base.clear();

    sk_base = vector<int>(param.n, 0L);
    sampler.get_binary_vector(sk_base);
    // sampler.get_ternary_hwt_vector(sk_ks, param.h);
    // cout << "Generation time of the secret key of the base scheme: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

void KeyGen::get_ksk_mk(KSKey_LWE_vec& ksk, const SKey_base_LWE& sk_base, const SKey_ks& sk_ks)
{
    ksk.resize(param.B_ksk_mk);
    int bound = param.B_ksk_mk >> 1;
    for (size_t bi = 0; bi < param.B_ksk_mk; bi++)
    {
        int a = (bi > bound) ? (bi - param.B_ksk_mk) : (bi);
        // reset key-switching key
        ksk[bi].A.clear();
        ksk[bi].b.clear();
        for (int i = 0; i < param.Nl_mk; i++)
        {
            vector<int> row(param.n,0L);
            ksk[bi].A.push_back(row);
        }
        ksk[bi].b = vector<int>(param.Nl_mk, 0L);

        vector<int> zg(param.Nl_mk, 0L);
        for (int i = 0; i < param.N; i++)
        {
            long sk_ks_i_long = long(a * sk_ks[i]);
            for (int j = 0; j < param.l_ksk_mk; j++)
            {
                long tmp_long = sk_ks_i_long * long(pow(param.B_ksk_mk, j));
                zg[param.l_ksk_mk * i + j] = parLWE.mod_q_base(tmp_long);
            }
        }

        // A as in the paper
        sampler.get_uniform_matrix(ksk[bi].A);

        // -A*s_base + e + zg
        normal_distribution<double> gaussian_sampler(0.0, Param::ksk_e_st_dev);
        for (int i = 0; i < param.Nl_mk; i++)
        {
            vector<int>& k_row = ksk[bi].A[i];
            for (int j = 0; j < param.n; j++)
                ksk[bi].b[i] -= k_row[j] * sk_base[j];
            ksk[bi].b[i] += (static_cast<int>(round(gaussian_sampler(rand_engine))) + zg[i]);
            parLWE.mod_q_base(ksk[bi].b[i]);
        }
    }

}

void KeyGen::get_bsk_mk_binary(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot)
{

    // clock_t start = clock();

    // index of a secret key coefficient of the base scheme
    int coef_counter = 0;

    // reset the input
    bsk.clear();

    // transform the secret key of the bootstrapping scheme to the DFT domain
    FFTPoly sk_boot_inv_fft(Param::N2p1);
    fftN.to_fft(sk_boot_inv_fft, sk_boot.sk_inv);

    // index of a secret key coefficient of the base scheme
    coef_counter = 0;
    // vector to keep the DFT transform of a random ternary vector
    FFTPoly g_fft(Param::N2p1);
    // vector to keep the DFT transform of a bootstrapping key part
    FFTPoly tmp_bsk_fft(Param::N2p1);
    // vector to keep a bootstrapping key part
    ModQPoly tmp_bsk(Param::N);
    vector<long> tmp_bsk_long;
    // precompute FFT transformed powers of decomposition bases
    vector<vector<FFTPoly>> B_bsk_pwr_poly;
    for (int iBase = 0; iBase < Param::B_bsk_size; iBase++)
    {
        double B_bsk_double = param.B_bsk[iBase];
        vector<FFTPoly> base_row;
        // FFT transform of (1,0,...,0)
        FFTPoly tmp_fft(Param::N2p1,complex<double>(1.0, 0.0));
        base_row.push_back(tmp_fft);
        for (int iPart = 1; iPart < param.l_bsk[iBase]; iPart++)
        {
            transform(tmp_fft.begin(), tmp_fft.end(), tmp_fft.begin(),
                        [B_bsk_double](complex<double> &z){ return z*B_bsk_double; });
            base_row.push_back(tmp_fft);
        }
        B_bsk_pwr_poly.push_back(base_row);
    }

    bsk.clear();
    bsk = vector<vector<NGSFFTctxt>>(Param::B_bsk_size);
    for (int i = 0; i < Param::B_bsk_size; i++)
        bsk[i] = vector<NGSFFTctxt>(param.bsk_partition[i], NGSFFTctxt(param.l_bsk[i], FFTPoly(Param::N2p1)));


    // loop over different decomposition bases
    for (int iBase = 0; iBase < Param::B_bsk_size; iBase++)
    {
        vector<NGSFFTctxt> base_row(param.bsk_partition[iBase], NGSFFTctxt(param.l_bsk[iBase], FFTPoly(Param::N2p1)));
        vector<FFTPoly>& B_bsk_pwr_poly_row = B_bsk_pwr_poly[iBase];
        for (int iCoef = coef_counter; iCoef < coef_counter+param.bsk_partition[iBase]; iCoef++)
        {
            NGSFFTctxt coef_row(param.l_bsk[iBase], FFTPoly(Param::N2p1));
            int sk_base_coef = sk_base[iCoef];
            // encrypt each bit using the NGS scheme
            for (int iPart = 0; iPart < param.l_bsk[iBase]; iPart++)
            {
                // sample random ternary vector
                ModQPoly g(Param::N,0L);
                sampler.get_ternary_hwt_vector(g, param.h);
                // FFT transform it
                fftN.to_fft(g_fft, g);
                // compute g * sk_boot^(-1)
                tmp_bsk_fft = g_fft* sk_boot_inv_fft;
                // compute g * sk_boot^(-1) + B^i * bit
                if (sk_base_coef == 1)
                    tmp_bsk_fft += B_bsk_pwr_poly_row[iPart];
                // inverse FFT of the above result
                fftN.from_fft(tmp_bsk_long, tmp_bsk_fft);
                // reduction modulo q_boot
                mod_q_boot(tmp_bsk, tmp_bsk_long);
                // FFT transform for further use
                fftN.to_fft(tmp_bsk_fft, tmp_bsk);

                coef_row[iPart] = tmp_bsk_fft;
            }
            base_row[iCoef-coef_counter] = coef_row;
        }
        bsk[iBase] = base_row;
        coef_counter += param.bsk_partition[iBase];
    }

    // cout << "Bootstrapping key generation: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}