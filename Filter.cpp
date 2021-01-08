#include "Filter.h"
#include "assert.h"
#include<math.h>

const double pi = 3.14159265;
const double e = 2.718281828459045;
void fftshift2D(cx_fmat* data, size_t xdim, size_t ydim);
void ifftshift2D(cx_fmat* data, size_t xdim, size_t ydim);
arma::fmat cx_abs(cx_fmat matrix);
arma::fmat real(arma::cx_fmat c_mat);
Filter::Filter()
{
}
void Filter::addData(fcube* img)
{
    // ´íÎó¼ì²é
    assert(rowHeightY > 1 && colWidthX > 1);
    this->img = img;
}
Mask::Mask(string maskType, int rowHeightY, int colWidthX)
{
    int cy = rowHeightY / 2;
    int cx = colWidthX / 2;
    maskArray = new cx_fmat(rowHeightY, colWidthX);
    if (maskType == "Threshold")
    {
        float threshold = 0.5;
        for (int r = 0; r < rowHeightY; r++) {
            for (int c = 0; c < colWidthX; c++) {
                if (sqrt((r - cy) * (r - cy) + (cx - c) * (cx - c)) <= (rowHeightY+colWidthX)/3) {
                    maskArray->at(r, c) = 1;
                }
                else {
                    maskArray->at(r, c) = 0;
                }
            }
        }
    }
    else if (maskType == "Gaussian")
    {
        double sig_square = 1;
        double base = sqrt(2 * pi) * sig_square;
        for (int r = 0; r < rowHeightY; r++) {
            for (int c = 0; c < colWidthX; c++) {
                double x = double(cx - c) / double(cx);
                double y = double(cy - r) / double(cy);
                x = x * x; y = y * y;
                double p = pow(e, -(x + y) / (2 * sig_square));
                maskArray->at(r, c) = p / base;
            }
        }
    }
}
Mask::~Mask()
{
    delete maskArray;
}
void Mask::domask(cx_fmat& shift)
{
    for (int r = 0; r < shift.n_rows; r++) {
        for (int c = 0; c < shift.n_cols; c++) {
            shift.at(r, c) = shift.at(r, c) * maskArray->at(r, c);
        }
    }
}
bool Filter::execute(string option, string maskType)
{

	if (option == "hpfilter") 
	{
        result = arma::fcube(img->n_rows, img->n_cols, img->n_slices);
        for (int s = 0; s < img->n_slices; s++) {
            cx_fmat img_fft = arma::fft2(img->slice(s));
            fftshift2D(&img_fft, img->n_rows, img->n_cols);
            Mask m(maskType, img->n_rows, img->n_cols);
            m.domask(img_fft);
            ifftshift2D(&img_fft, img->n_rows, img->n_cols);
            cx_fmat img_ifft = ifft2(img_fft);
            result.slice(s) = cx_abs(img_ifft);
            return true;
        }
	}
	else if (option == "lpfilter") 
	{
		return false;
	}
	else 
	{
		return false;
	}
}


Filter::~Filter()
{
	//delete img;
}

void fftshift2D(cx_fmat* data, size_t xdim, size_t ydim)
{
    size_t xshift = xdim / 2;
    size_t yshift = ydim / 2;
    if ((xdim * ydim) % 2 != 0) {
        std::vector<cx_float> out;
        out.resize(xdim * ydim);
        for (size_t x = 0; x < xdim; x++) {
            size_t outX = (x + xshift) % xdim;
            for (size_t y = 0; y < ydim; y++) {
                size_t outY = (y + yshift) % ydim;
                out[outX + xdim * outY] = data->at(x, y);
            }
        }
        copy(out.begin(), out.end(), &data[0]);
    }
    else {
        for (size_t x = 0; x < xdim; x++) {
            size_t outX = (x + xshift) % xdim;
            for (size_t y = 0; y < yshift; y++) {
                size_t outY = (y + yshift) % ydim;
                swap(data->at(outX, outY), data->at(x, y));
                //swap(data[outX + xdim * outY], data[x + xdim * y]);
            }
        }
    }
}


void ifftshift2D(cx_fmat* data, size_t xdim, size_t ydim)
{
    size_t xshift = xdim / 2;
    if (xdim % 2 != 0) {
        xshift++;
    }
    size_t yshift = ydim / 2;
    if (ydim % 2 != 0) {
        yshift++;
    }
    if ((xdim * ydim) % 2 != 0) {
        std::vector<cx_float > out;
        out.resize(xdim * ydim);
        for (size_t x = 0; x < xdim; x++) {
            size_t outX = (x + xshift) % xdim;
            for (size_t y = 0; y < ydim; y++) {
                size_t outY = (y + yshift) % ydim;
                out[outX + xdim * outY] = data->at(x, y);
            }
        }
        copy(out.begin(), out.end(), &data[0]);
    }
    else {
        for (size_t x = 0; x < xdim; x++) {
            size_t outX = (x + xshift) % xdim;
            for (size_t y = 0; y < yshift; y++) {
                size_t outY = (y + yshift) % ydim;
                swap(data->at(outX, outY), data->at(x, y));
                //swap(data[outX + xdim * outY], data[x + xdim * y]);
            }
        }
    }
}
arma::fmat real(arma::cx_fmat c_mat)
{
    arma::fmat mat(c_mat.n_rows, c_mat.n_cols);

    for (int r = 0; r < c_mat.n_rows; r++) {
        for (int c = 0; c < c_mat.n_cols; c++) {
            mat.at(r, c) = c_mat.at(r, c).real();
        }
    }
    return mat;
}
arma::fmat cx_abs(cx_fmat matrix)
{
    arma::fmat result(matrix.n_rows, matrix.n_cols);
    for (int r = 0; r < matrix.n_rows; r++) {
        for (int c = 0; c < matrix.n_cols; c++) {
            float real_val = matrix.at(r, c).real();
            float img_val = matrix.at(r, c).imag();
            result.at(r, c) = sqrt(real_val * real_val + img_val * img_val);
        }
    }
    return result;
}