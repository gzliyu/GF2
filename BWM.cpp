#include "BWM.h"
arma::vec cx_real(arma::cx_vec c_v)
{
    arma::vec v(c_v.n_rows);
    for (int i = 0; i < c_v.n_rows; i++) {
        v.at(i) = c_v.at(i).real();
    }
    return v;
}
arma::fmat cx_real(arma::cx_fmat c_mat)
{
    arma::fmat mat(c_mat.n_rows, c_mat.n_cols);

    for (int r = 0; r < c_mat.n_rows; r++) {
        for (int c = 0; c < c_mat.n_cols; c++) {
            mat.at(r, c) = c_mat.at(r, c).real();
        }
    }
    return mat;
}
int Rand(int n)
{
    return rand() % n;
}
arma::cube BWM::encode(arma::cube& x, arma::mat& wm)
{
    arma::cube x_result(x.n_rows, x.n_cols, x.n_slices);
    arma::cx_cube x_fft(x.n_rows, x.n_cols, x.n_slices);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    for (int s = 0; s < x.n_slices; s++) {
        arma::cx_mat img_mat_fft = fft2(x.slice(s));
        x_fft.slice(s) = img_mat_fft;

        std::vector<int>rlist;
        for (int r = 0; r<int(x.n_rows / 2); r++)
            rlist.push_back(r);
        std::vector<int>clist;
        for (int c = 0; c < x.n_cols; c++)
            clist.push_back(c);
        std::random_shuffle(rlist.begin(), rlist.end(), std::pointer_to_unary_function<int, int>(Rand));
        std::random_shuffle(clist.begin(), clist.end(), std::pointer_to_unary_function<int, int>(Rand));
        arma::mat tmp(x.n_rows, x.n_cols);
        for (int r = 0; r<int(x.n_rows / 2); r++) {
            for (int c = 0; c < x.n_cols; c++) {
                if (rlist[r] < wm.n_rows && clist[c] < wm.n_cols) {
                    tmp.at(r, c) = wm.at(rlist[r], clist[c]);
                    tmp.at(x.n_rows - 1 - r, x.n_cols - 1 - c) = tmp.at(r, c);
                }
            }
        }
        arma::cx_mat result(x.n_rows, x.n_cols);
        float alpha = 2000;
        result = alpha * tmp + x_fft.slice(s);
        arma::cx_mat result_mat_ifft = ifft2(result);
        x_result.slice(s) = cx_real(result_mat_ifft);
    }
    return x_result;
}

arma::cube BWM::encode(arma::fcube& x, arma::mat& wm)
{
    arma::cube x_conv = arma::conv_to<arma::cube>::from(x);
    return encode(x_conv, wm);
}



arma::fcube BWM::encode(arma::fcube& x, arma::fmat& wm)
{
    arma::fcube x_result(x.n_rows, x.n_cols, x.n_slices);
    arma::cx_fcube x_fft(x.n_rows, x.n_cols, x.n_slices);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    for (int s = 0; s < x.n_slices; s++) {
        arma::cx_fmat img_mat_fft = fft2(x.slice(s));
        x_fft.slice(s) = img_mat_fft;
        std::vector<int>rlist;
        for (int r = 0; r<int(x.n_rows / 2); r++)
            rlist.push_back(r);
        std::vector<int>clist;
        for (int c = 0; c < x.n_cols; c++)
            clist.push_back(c);
        std::random_shuffle(rlist.begin(), rlist.end(), std::pointer_to_unary_function<int, int>(Rand));
        std::random_shuffle(clist.begin(), clist.end(), std::pointer_to_unary_function<int, int>(Rand));
        arma::fmat tmp(x.n_rows, x.n_cols);
        for (int r = 0; r<int(x.n_rows / 2); r++) {
            for (int c = 0; c < x.n_cols; c++) {
                if (rlist[r] < wm.n_rows && clist[c] < wm.n_cols) {
                    tmp.at(r, c) = wm.at(rlist[r], clist[c]);
                    tmp.at(x.n_rows - 1 - r, x.n_cols - 1 - c) = tmp.at(r, c);
                }
            }
        }
        arma::cx_fmat result(x.n_rows, x.n_cols);
        float alpha = 2000;
        result = alpha * tmp + x_fft.slice(s);
        arma::cx_fmat result_mat_ifft = ifft2(result);
        x_result.slice(s) = cx_real(result_mat_ifft);
    }
    printf("complete!\n");
    return x_result;
}
