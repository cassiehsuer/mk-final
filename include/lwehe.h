#ifndef LWEHE
#define LWEHE

#include "params.h"
#include "keygen.h"

class Ctxt_LWE
{
    public:
    std::vector<vector<int>> a;
    int b;

    Ctxt_LWE()
    {
        a.clear();
        a.resize(1, vector<int>(parLWE.n, 0L));
    }
    Ctxt_LWE(const Ctxt_LWE& ct);

    Ctxt_LWE& operator=(const Ctxt_LWE& ct);

    Ctxt_LWE operator +(const Ctxt_LWE& ct) const;
    Ctxt_LWE operator -(const Ctxt_LWE& ct) const;
};

Ctxt_LWE operator -(const int c, const Ctxt_LWE& ct);

/**
 * Switches a given ciphertext to a given modulus.
 * @param[in,out] ct ciphertext
 * @param[in] old_q old modulus
 * @param[in] new_q old modulus
 */
inline void modulo_switch_lwe(Ctxt_LWE& ct, int old_q, int new_q)
{
    for (size_t i = 0; i < ct.a.size(); i++)
    {
        std::vector<int>& a = ct.a[i];
        for (size_t j = 0; j < a.size(); j++)
            a[j] = int((a[j]*new_q)/old_q);
    }
    ct.b = int((ct.b*new_q)/old_q);
}

/**
 * Switches a given polynomial from q_base to modulus 2*N.
 * @param[in,out] poly polynomial
 */
inline void modulo_switch_to_boot(Ctxt_LWE& poly)
{
    modulo_switch_lwe(poly, parLWE.q_base, Param::N2);
}

/**
 * Switches a given polynomial from q_boot to q_base.
 * @param[in,out] poly polynomial
 */
inline void modulo_switch_to_base_lwe(ModQPoly& poly)
{
    modulo_switch(poly, q_boot, parLWE.q_base);
}

inline void modulo_switch_to_base_lwe_vec(Ctxt_LWE& ct, int old_q, int new_q)
{
    double ratio = double(new_q)/double(old_q);
    for (size_t i = 0; i < ct.a.size(); i++)
    {
        std::vector<int>& a = ct.a[i];
        for (size_t j = 0; j < a.size(); j++)
            a[j] = int(round(a[j]*ratio));
    }
    ct.b = int(round(ct.b*ratio));
}

/**
 * Computes the external product of a given polynomial ciphertext
 * with an NGS ciphertext in the FFT form
 * @param[in,out] poly polynomial ciphertext
 * @param[in] poly_vector NGS ciphertext
 * @param[in] b decomposition base, power of 2
 * @param[in] shift bit shift to divide by b
 * @param[in] l decomposition length
 */
void external_product(std::vector<long>& res, const std::vector<int>& poly, const std::vector<FFTPoly>& poly_vector, int b, int shift, int l);

class SchemeLWE
{


    public:

        SchemeLWE()
        {
        }


        /**
         * Encrypts a bit using LWE.
         * @param[out] ct ciphertext encrypting the input bit
         * @param[in] m bit to encrypt
         */

        void encrypt_mk(Ctxt_LWE& ct, int m, vector<SKey_base_LWE> sk_base_all) const;

        /**
         * Decrypts a bit using LWE.
         * @param[out] ct ciphertext encrypting a bit
         * @return b bit
         */
        int decrypt(const Ctxt_LWE& ct, vector<SKey_base_LWE> sk_base_all, int t) const;

        /**
         * Performs key switching of a given ciphertext
         * to LWE
         * @param[out] ct LWE ciphertext (vector of dimension n)
         * @param[in] poly polynomial ciphertext (vector of dimension N)
         */

        void key_switch_mk(Ctxt_LWE& ct_out, const Ctxt_LWE& ct_in, vector<KSKey_LWE_vec>& ksk_all) const;


        void mk_bootstrap_binary(Ctxt_LWE& ct, vector<BSKey_LWE> bsk_all, vector<HybridKey>& hybridkey_all, vector<KSKey_LWE_vec>& ksk_all) const;


        /**
         * Computes the NAND gate of two given ciphertexts ct1 and ct2
         * @param[out] ct_res encryptions of the outuput of the NAND gate
         * @param[in] ct_1 encryption of the first input bit
         * @param[in] ct_2 encryption of the second input bit
         */
        void nand_gate(Ctxt_LWE& ct_res, const Ctxt_LWE& ct1, const Ctxt_LWE& ct2) const;


};





#endif