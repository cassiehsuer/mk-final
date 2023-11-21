#include "lwehe.h"
#include "sampler.h"
#include "fft.h"

#include <iostream>
#include <time.h>
#include <chrono>
#include <vector>
#include <cassert>


using namespace std;

Ctxt_LWE::Ctxt_LWE(const Ctxt_LWE& ct)
{
    a = ct.a;
    b = ct.b;
}

Ctxt_LWE& Ctxt_LWE::operator=(const Ctxt_LWE& ct)
{
    a = ct.a;
    b = ct.b;
    return *this;
}

Ctxt_LWE Ctxt_LWE::operator +(const Ctxt_LWE& ct) const
{
    Ctxt_LWE res;
    res.a = vector<vector<int>>(ct.a.size());
    for (size_t i = 0; i < ct.a.size(); i++)
    {
        res.a[i].resize(parLWE.n);
        for (size_t j = 0; j < parLWE.n; j++)
            res.a[i][j] = parLWE.mod_q_base(a[i][j] + ct.a[i][j]);
    }
    res.b = parLWE.mod_q_base(b + ct.b);
    return res;
}

Ctxt_LWE Ctxt_LWE::operator -(const Ctxt_LWE& ct) const
{
    Ctxt_LWE res;
    res.a = vector<vector<int>>(ct.a.size());
    for (size_t i = 0; i < ct.a.size(); i++)
    {
        res.a[i].resize(parLWE.n);
        for (size_t j = 0; j < parLWE.n; j++)
            res.a[i][j] = parLWE.mod_q_base(a[i][j] - ct.a[i][j]);
    }
    res.b = parLWE.mod_q_base(b - ct.b);
    return res;
}

Ctxt_LWE operator -(const int c, const Ctxt_LWE& ct)
{
    Ctxt_LWE res;
    res.a = vector<vector<int>>(ct.a.size());
    for (size_t i = 0; i < ct.a.size(); i++)
    {
        res.a[i].resize(parLWE.n);
        for (size_t j = 0; j < parLWE.n; j++)
            res.a[i][j] = parLWE.mod_q_base(-ct.a[i][j]);
    }
    res.b = parLWE.mod_q_base(c-ct.b);
    return res;
}

void SchemeLWE::encrypt_mk(Ctxt_LWE& ct, int m, vector<SKey_base_LWE> sk_base_all) const
{
    clock_t start = clock();

    int n = parLWE.n;

    normal_distribution<double> gaussian_sampler(0.0, Param::e_st_dev);
    int b = parLWE.delta_base*m + static_cast<int>(round(gaussian_sampler(rand_engine)));
    Sampler s(parLWE);
    ct.a.clear();
    ct.a.resize(sk_base_all.size(), vector<int>(parLWE.n, 0L));

    for (int j = 0; j < sk_base_all.size(); j++){
        s.get_uniform_vector(ct.a[j]);
        for (int i = 0; i < n; i++)
        {
            b -= sk_base_all[j][i] * ct.a[j][i];
        }
    }
    parLWE.mod_q_base(b);
    ct.b = b;

    // cout << "Encryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}

int SchemeLWE::decrypt(const Ctxt_LWE& ct, vector<SKey_base_LWE> sk_base_all, int t) const
{
    // clock_t start = clock();
    int output = ct.b;
    int n = parLWE.n;
    for (int i = 0; i < sk_base_all.size(); i++)
    {
        for (int j = 0; j < n; j++)
        {
            output += ct.a[i][j] * sk_base_all[i][j];
        }
    }
    output = parLWE.mod_q_base(output);
    output = int(round(double(output*t)/double(parLWE.q_base)));
    //cout << "Decryption: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
    return output;
}

void external_product(vector<long>& res, const vector<int>& poly, const vector<FFTPoly>& poly_vector, int b, int shift, int l)
{
    int N = Param::N;
    int N2p1 = Param::N2p1;

    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    vector<int> poly_decomp(N);

    for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    FFTPoly res_fft(N2p1);
    FFTPoly tmp_fft(N2p1);
    int mask = b-1;
    int bound = b >> 1;
    int digit, sgn;
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int& abs_val = poly_abs[i];
            digit = abs_val & mask;
            if (digit > bound)
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? (digit - b): (b - digit);
                abs_val >>= shift;
                ++abs_val;
            }
            else
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? digit: -digit;
                abs_val >>= shift;
            }
        }
        fftN.to_fft(tmp_fft, poly_decomp);
        tmp_fft *= poly_vector[j];
        res_fft += tmp_fft;
    }
    fftN.from_fft(res, res_fft);
}

void SchemeLWE::key_switch_mk(Ctxt_LWE& ct_out, const Ctxt_LWE& ct_in, vector<KSKey_LWE_vec>& ksk_all) const{
    int k = ksk_all.size();
    int N = Param::N;
    int n = parLWE.n;
    int Nl_mk = parLWE.Nl_mk;
    int B_ksk_mk = Param::B_ksk_mk;
    int l_ksk_mk = parLWE.l_ksk_mk;
    long b_i_out;
    int bound = B_ksk_mk >> 1;
    int digit;
    int jl;
    int tmp;
    int sgn;
    int bi;
    ct_out.b = ct_in.b;
    for (int i = 0; i < k; i++)
    {
        vector<int> a_decomp(Nl_mk);
        vector<int> a_i_sign(N);
        vector<int> a_i_abs(N);
        vector<long> a_i_out(n, 0L);
        b_i_out = 0L;
        jl = 0;
        const vector<int>& a_i = ct_in.a[i];
        for (int j = 0; j < N; ++j)
        {
            const int& aij = a_i[j];
            a_i_abs[j] = abs(aij);
            a_i_sign[j] = (aij < 0)? -1 : 1;
        }
        for (int j = 0; j < N; ++j)
        {
            tmp = a_i_abs[j];
            sgn = a_i_sign[j];
            for (int t = 0; t < l_ksk_mk; ++t)
            {
                digit = tmp % B_ksk_mk;
                if (digit > bound)
                {
                    a_decomp[jl+t] = (sgn == 1) ? (digit - B_ksk_mk) : (B_ksk_mk - digit);
                    tmp /= B_ksk_mk;
                    ++tmp;
                }
                else
                {
                    a_decomp[jl+t] = (sgn == 1) ? digit: - digit;
                    tmp /= B_ksk_mk;
                }
                bi = (a_decomp[jl+t] < 0) ? (a_decomp[jl+t] + B_ksk_mk) : a_decomp[jl+t];
                const vector<int>& ksk_row = ksk_all[i][bi].A[jl+t];
                for (int f = 0; f < n; f++)
                {
                    a_i_out[f] += long(ksk_row[f]);
                }
                const vector<int>& ksk_b = ksk_all[i][bi].b;
                b_i_out += long(ksk_b[jl+t]);
            }
            jl += l_ksk_mk;
        }
        parLWE.mod_q_base(ct_out.a[i], a_i_out);
        ct_out.b += parLWE.mod_q_base(b_i_out);
    }
}

void SchemeLWE::mk_bootstrap_binary(Ctxt_LWE& ct_n, vector<BSKey_LWE> bsk_all, vector<HybridKey>& hybridkey_all, vector<KSKey_LWE_vec>& ksk_all) const {

    int N = Param::N;
    int N2 = Param::N2;
    int N2p1 = Param::N2p1;
    int B_bsk_size = Param::B_bsk_size;
    int half_delta_boot = Param::half_delta_boot;
    int k = bsk_all.size();

    // switch to modulus 2*N
    modulo_switch_to_boot(ct_n);

    // initialize accumulator by Q/8*X^b*X^{N/2}*(1+X+...+X^{N-1})
    vector<int> acc(N, half_delta_boot);
    int b_sign = 1;
    int b_pow = (N/2 + ct_n.b) % N2;
    if (b_pow < 0)
        b_pow += N2;
    if (b_pow >= N)
    {
        b_pow -= N;
        b_sign = -1;
    }
    for (int i = 0; i < b_pow; ++i)
        acc[i] = (b_sign == 1) ? -acc[i]: acc[i];
    for (int i = b_pow; i < N; ++i)
        acc[i] = (b_sign == 1) ? acc[i]: -acc[i];

    //accumulator loop
    int coef, coef_sign, B, shift, l;
    double Bd;
    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);

    for (size_t iBsk = 0; iBsk < k; ++iBsk)
    {
        int coef_counter = 0;
        for (int iBase = 0; iBase < B_bsk_size; ++iBase)
        {
            B = parLWE.B_bsk[iBase];
            Bd = double(B);
            shift = parLWE.shift_bsk[iBase];
            l = parLWE.l_bsk[iBase];
            const vector<NGSFFTctxt>& bk_coef_row = bsk_all[iBsk][iBase];
            for (int iCoef = 0; iCoef < parLWE.bsk_partition[iBase]; ++iCoef)
            {
                coef = ct_n.a[iBsk][iCoef+coef_counter];
                if (coef == 0) continue;
                coef_sign = 1;
                if (coef < 0) coef += N2;
                if (coef >= N)
                {
                    coef -= N;
                    coef_sign = -1;
                }

                // acc * (X^coef - 1)
                if (coef_sign == 1)
                {
                    for (int i = 0; i<coef; ++i)
                        tmp_poly[i] = mod_q_boot(-acc[i-coef+N] - acc[i]);
                    for (int i = coef; i < N; ++i)
                        tmp_poly[i] = mod_q_boot(acc[i-coef] - acc[i]);
                }
                else
                {
                    for (int i = 0; i<coef; ++i)
                        tmp_poly[i] = mod_q_boot(acc[i-coef+N] - acc[i]);
                    for (int i = coef; i < N; ++i)
                        tmp_poly[i] = mod_q_boot(-acc[i-coef] - acc[i]);
                }

                // acc * (X^coef - 1) x bk[i]
                external_product(tmp_poly_long, tmp_poly, bk_coef_row[iCoef], B, shift, l);
                mod_q_boot(tmp_poly, tmp_poly_long);
                // acc * (X^coef - 1) x bk[i] + acc
                for (int i = 0; i<N; ++i)
                    acc[i] += tmp_poly[i];
            }

            coef_counter += parLWE.bsk_partition[iBase];
        }
    }
    // add floor(q_boot/(2*t)) to all coefficients of the accumulator
    for (auto it = acc.begin(); it < acc.end(); ++it)
        *it += half_delta_boot;

    //mod q_boot of the accumulator
    mod_q_boot(acc);

    vector<ModQPoly> c_hybrid(k, ModQPoly(N,0));
    c_hybrid.push_back(acc);

    vector<vector<long>> u_long(k + 1, vector<long>(N));
    vector<vector<long>> v_long(k + 1, vector<long>(N));
    vector<vector<long>> w_0_long(k + 1, vector<long>(N));
    vector<vector<long>> w_1_long(k + 1, vector<long>(N));
    vector<vector<int>> u_int(k + 1, vector<int>(N));
    vector<vector<int>> v_int(k + 1, vector<int>(N));
    vector<vector<int>> w_0_int(k + 1, vector<int>(N));
    vector<vector<int>> w_1_int(k + 1, vector<int>(N));
    int B_hybrid = parLWE.B_hybrid;
    int shift_hybrid = parLWE.shift_hybrid;
    int l_hybrid = parLWE.l_hybrid;

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j <= k; j++)
        {
            external_product(u_long[j], c_hybrid[j], hybridkey_all[i].d_2, B_hybrid, shift_hybrid, l_hybrid);
            mod_q_boot(u_int[j], u_long[j]);

            external_product(v_long[j], c_hybrid[j], hybridkey_all[j].b, B_hybrid, shift_hybrid, l_hybrid);
            mod_q_boot(v_int[j], v_long[j]);

            external_product(w_0_long[j], v_int[j], hybridkey_all[i].d_0, B_hybrid, shift_hybrid, l_hybrid);
            mod_q_boot(w_0_int[j], w_0_long[j]);

            external_product(w_1_long[j], v_int[j], hybridkey_all[i].d_1, B_hybrid, shift_hybrid, l_hybrid);
            mod_q_boot(w_1_int[j], w_1_long[j]);
        }
        c_hybrid.clear();
        c_hybrid.resize(k + 1, ModQPoly(Param::N,0));
        for (int j = 0; j <= k; j++)
        {
            for (size_t t = 0; t < N; t++)
            {
                c_hybrid[j][t] += u_int[j][t];
                c_hybrid[k][t] += w_0_int[j][t];
                c_hybrid[i][t] += w_1_int[j][t];
            }
            mod_q_boot(c_hybrid[j]);
        }
    }

    Ctxt_LWE ct_NN;
    ct_NN.b = c_hybrid[k][0];
    ct_NN.a.clear();
    ct_NN.a.resize(k, vector<int>(N, 0L));
    for (int i = 0; i < k; i++)
    {
        ct_NN.a[i][0] = c_hybrid[i][0];
        for (size_t j = 1; j < N; j++)
        {
            ct_NN.a[i][j] = -(c_hybrid[i][N-j]);
        }
    }

    modulo_switch_to_base_lwe_vec(ct_NN, q_boot, parLWE.q_base);

    // clock_t start_ks = clock();
    key_switch_mk(ct_n, ct_NN, ksk_all);
    // cout<<"KS time:"<<float(clock()-start_ks)/CLOCKS_PER_SEC;

}


void SchemeLWE::nand_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const
{
    //clock_t start = clock();
    ct_res = parLWE.nand_const - (ct1 + ct2);
    //cout << "NAND: " << float(clock()-start)/CLOCKS_PER_SEC << endl;
}
