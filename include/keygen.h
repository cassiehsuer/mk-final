#ifndef KEYGEN
#define KEYGEN

#include <NTL/mat_ZZ.h>
#include <vector>
#include <fftw3.h>
#include <complex>
#include "params.h"
#include "sampler.h"

// secret key of the bootstrapping scheme
typedef struct {
    ModQPoly sk;
    ModQPoly sk_inv;
} SKey_boot;

// Hybrid product public key
typedef struct {
    NGSFFTctxt b;
    NGSFFTctxt d_1;
    NGSFFTctxt d_0;
    NGSFFTctxt d_2;
} HybridKey;

// secret key of the LWE base scheme
typedef std::vector<int> SKey_base_LWE;

typedef std::vector<int> SKey_ks;

// secret key of the hybrid product
typedef ModQPoly SKey_hybrid;

/**
 * Bootstrapping key.
 * It consists of several sets of keys corresponding to different
 * decomposition bases of the bootstrapping key B_bsk.
 * The i-th set contains vectors with l_bsk[i] complex vectors.
 * These complex vectors are an encryption of some bit of the secret key
 * of the base scheme in the NGS form.
 */
typedef std::vector<std::vector<NGSFFTctxt>> BSKey_LWE;

// key-switching key
typedef struct{
    ModQMatrix A;
    std::vector<int> b;
}  KSKey_LWE;

typedef std::vector<KSKey_LWE> KSKey_LWE_vec;

class KeyGen
{
    Param param;
    Sampler sampler;

    public:

        KeyGen(Param _param): param(_param), sampler(_param)
        {}

        // generate common crs a for hybrid product
        void get_crs_a(NGSFFTctxt& crs_a, NGSFFTctxt& crs_a_neg);

        // generate the sk of hybrid product
        void get_sk_hybrid(SKey_hybrid& sk_hybrid);

        // generate the public key of hybrid product
        void get_hybridkey(HybridKey& hybridkey, SKey_hybrid& sk_hybrid, SKey_boot& sk_boot, NGSFFTctxt& crs_a);

        /**
         * Generate a secret key of the bootstrapping scheme.
         * @param[out] sk_boot secret key of the bootstrapping scheme.
         */
        void get_sk_boot(SKey_boot& sk_boot);

        /**
         * Generate a secret key of the base scheme.
         * @param[out] sk_base secret key of the base scheme.
         */
        void get_sk_base(SKey_base_LWE& sk_base);

        void get_sk_ks(SKey_ks& sk_ks);

        void get_ksk_mk(KSKey_LWE_vec& ksk, const SKey_base_LWE& sk_base, const SKey_ks& sk_ks);

        void get_bsk_mk_binary(BSKey_LWE& bsk, const SKey_base_LWE& sk_base, const SKey_boot& sk_boot);

};

#endif