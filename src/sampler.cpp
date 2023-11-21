#include "sampler.h"
#include "params.h"
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>
#include <cassert>

using namespace std;

void Sampler::get_ternary_hwt_vector(vector<int>& vec, int h)
{
    int idx = 0;
    uniform_int_distribution<int> vec_size_sampler(0, vec.size() - 1);
	while(idx < h) {
		int i = vec_size_sampler(rand_engine);
		if(vec[i] == 0) {
			vec[i] = (binary_sampler(rand_engine) == 0) ? 1 : -1;
			idx ++;
		}
	}
}

void Sampler::get_ternary_vector(vector<int>& vec)
{

    for(int i=0; i<vec.size(); i++)
        vec[i] = ternary_sampler(rand_engine);
}

void Sampler::get_ternary_matrix(vector<vector<int>>& mat)
{

    for(int i=0; i<mat.size(); i++)
    {
        vector<int>& row = mat[i];
        get_ternary_vector(row);
    }
}

void Sampler::get_binary_vector(vector<int>& vec)
{
    for(int i=0; i<vec.size(); i++)
        vec[i] = binary_sampler(rand_engine);
}

void Sampler::get_uniform_vector(vector<int>& vec)
{
    for(int i=0; i<vec.size(); i++)
        vec[i] = mod_q_base_sampler(rand_engine);
}

void Sampler::get_uniform_vector_boot(vector<int>& vec)
{
    for(int i=0; i<vec.size(); i++)
        vec[i] = mod_q_boot_sampler(rand_engine);
}

void Sampler::get_uniform_vector_with_neg(vector<int>& vec, vector<int>& vec_neg)
{
    for(int i=0; i<vec.size(); i++){
        vec[i] = mod_q_boot_sampler(rand_engine);
        vec_neg[i] = - vec[i];
    }
}

void Sampler::get_uniform_matrix(vector<vector<int>>& mat)
{
    for(int i=0; i<mat.size(); i++)
    {
        vector<int>& row = mat[i];
        get_uniform_vector(row);
    }
}

void Sampler::get_uniform_matrix_boot(vector<vector<int>>& mat)
{
    for(int i=0; i<mat.size(); i++)
    {
        vector<int>& row = mat[i];
        get_uniform_vector_boot(row);
    }
}

void Sampler::get_uniform_matrix_with_neg(vector<vector<int>>& mat, vector<vector<int>>& mat_neg)
{
    for(int i=0; i<mat.size(); i++)
    {
        vector<int>& row = mat[i];
        vector<int>& row_neg = mat_neg[i];
        get_uniform_vector_with_neg(row, row_neg);
    }
}

void Sampler::get_gaussian_vector(vector<int>& vec, double st_dev)
{
    normal_distribution<double> gaussian_sampler(0.0, st_dev);
    for(size_t i=0; i<vec.size(); i++)
        vec[i] = static_cast<int>(round(gaussian_sampler(rand_engine)));
}

void Sampler::get_gaussian_matrix(vector<vector<int>>& mat, double st_dev)
{
    for(size_t i=0; i<mat.size(); i++)
    {
        vector<int>& row = mat[i];
        get_gaussian_vector(row, st_dev);
    }
}

void Sampler::get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale, int shift)
{
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    //element of Z_(q_boot)
    ZZ_p coef;
    coef.init(ZZ(q_boot));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling
    while (true)
    {
        //create the polynomial with the coefficient vector of the desired form

        vector<int> tmp_vec(param.N, 0L);
        get_ternary_hwt_vector(tmp_vec, param.h);
        SetCoeff(poly, 0, tmp_vec[0]*scale + shift);
        for (size_t i = 1; i < vec.size(); i++)
        {
            coef = tmp_vec[i]*scale;
            SetCoeff(poly, i, coef);
        }
        //test invertibility
        try
        {
            InvMod(inv_poly, poly, Param::get_def_poly());
            break;
        }
        catch(...)
        {
            cout << "Polynomial " << poly << " isn't a unit" << endl;
            continue;
        }
    }
    //cout << "Poly: " << poly << endl;
    //cout << "Poly inverse: " << inv_poly << endl;
    //extract the coefficient vector of poly
    int tmp_coef;
    for (int i = 0; i <= deg(poly); i++)
    {
        tmp_coef = conv<long>(poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (int i = 0; i <= deg(inv_poly); i++)
    {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }

    //cout << "Vector:" << vec << endl;
    //cout << "Inverse vector:" << vec_inv << endl;
}

void Sampler::get_invertible_matrix(vector<vector<int>>& mat, vector<vector<int>>& mat_inv, int scale, int shift)
{
    //check that the input matrices are squares
    assert(mat[0].size() == mat.size());
    assert(mat_inv[0].size() == mat_inv.size());
    //check that both input matrices have the same dimension
    assert(mat.size() == mat_inv.size());

    //number of rows of the input matrix
    int dim = mat.size();

    //element of Z_(q_boot)
    ZZ_p coef;
    coef.init(ZZ(param.q_base));

    //candidate matrix
    mat_ZZ_p tmp_mat(INIT_SIZE, dim, dim);

    //candidate inverse matrix
    mat_ZZ_p tmp_mat_inv(INIT_SIZE, dim, dim);

    //sampling and testing
    while (true)
    {
        //sampling
        for (int i = 0; i < dim; i++)
        {
            Vec<ZZ_p>& row = tmp_mat[i];
            for (int j = 0; j < dim; j++)
            {
                coef = ternary_sampler(rand_engine)*scale;
                if (i==j)
                    coef += ZZ_p(shift);
                row[j] = coef;
            }
        }
        //test invertibility
        try
        {
            inv(tmp_mat_inv, tmp_mat);
            break;
        }
        catch(...)
        {
            cout << "Matrix " << tmp_mat << " is singular" << endl;
            continue;
        }
    }
    //lift mod q representation to integers
    int tmp_coef;
    for (int i = 0; i < dim; i++)
    {
        Vec<ZZ_p>& tmp_row = tmp_mat[i];
        Vec<ZZ_p>& tmp_row_inv = tmp_mat_inv[i];
        vector<int>& row = mat[i];
        vector<int>& row_inv = mat_inv[i];
        for (int j = 0; j < dim; j++)
        {
            tmp_coef = conv<long>(tmp_row[j]);
            if (tmp_coef > param.half_q_base)
                tmp_coef -= param.q_base;
            row[j] = tmp_coef;

            tmp_coef = conv<long>(tmp_row_inv[j]);
            if (tmp_coef > param.half_q_base)
                tmp_coef -= param.q_base;
            row_inv[j] = tmp_coef;
        }
    }
}