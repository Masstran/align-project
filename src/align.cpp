#include "align.h"
#include <string>
#include <limits>
#include <cmath>

using std::string;
using std::cout;
using std::endl;

Image mirror(const Image src_image, const uint radius)
{
    Image tmp(src_image.n_rows + 2*radius, src_image.n_cols + 2*radius);
    for (uint i = 0; i < src_image.n_rows + 2*radius; i++)
    {
        for (uint j = 0; j < src_image.n_cols + 2*radius; j++)
        {
            int tmpi = (i < radius) ? radius - i - 1 : (i >= src_image.n_rows + radius) ? 2*src_image.n_rows + radius - i - 1 : i - radius;
            int tmpj = (j < radius) ? radius - j - 1 : (j >= src_image.n_cols + radius) ? 2*src_image.n_cols + radius - j - 1 : j - radius;
            tmpi = (tmpi < 0) ? 0 : (static_cast<uint>(tmpi) >= src_image.n_rows) ? src_image.n_rows - 1 : tmpi;
            tmpj = (tmpj < 0) ? 0 : (static_cast<uint>(tmpj) >= src_image.n_cols) ? src_image.n_cols - 1 : tmpj;
            //uint tmpj = (j == 0) ? 0 : (j == out.n_cols + 1) ? out.n_cols - 1 : j - 1;
            //std::cout << i << ' ' << tmpi << ' ' << j << ' ' << tmpj << '\n';
            tmp(i, j) = src_image(tmpi, tmpj);
        }
    }
    return tmp;
}

double MSE(const Image I1, const Image I2, const int i_off, const int j_off)
{
    double res = 0;
    int sz = 0;
    uint i1, i2;
    const uint r = I1.n_rows/10, c = I1.n_cols/10;
    //const uint r = 0, c = 0;
    for (int i = r; static_cast<unsigned>(i) < I1.n_rows - r; i++)
    {
        for (int j = c; static_cast<unsigned>(j) < I1.n_cols - c; j++)
        {
            if ((i + i_off >= 0) && (static_cast<unsigned>(i + i_off) < I1.n_rows) &&
                (j + j_off >= 0) && (static_cast<unsigned>(j + j_off) < I1.n_cols))
            {
                sz++;
                i1 = std::get<0>(I1(i + i_off, j + j_off));
                i2 = std::get<0>(I2(i, j));
                res += (i1 - i2)*(i1 - i2);
            }
        }
    }
    if (sz == 0)
        return std::numeric_limits<double>::infinity();
    res /= static_cast<double>(sz);
    return res;
}

double cross_corr(const Image I1, const Image I2, const int i_off, const int j_off)
{
    double res = 0;
    //int sz = 0;
    uint i1, i2;
    const uint r = I1.n_rows/10, c = I1.n_cols/10;
    for (int i = r; static_cast<unsigned>(i) < I1.n_rows - r; i++)
    {
        for (int j = c; static_cast<unsigned>(j) < I1.n_cols - c; j++)
        {
            if ((i + i_off >= 0) && (static_cast<unsigned>(i + i_off) < I1.n_rows) &&
                (j + j_off >= 0) && (static_cast<unsigned>(j + j_off) < I1.n_cols))
            {
                //sz++;
                i1 = std::get<0>(I1(i + i_off, j + j_off));
                i2 = std::get<0>(I2(i, j));
                res += i1*i2;
            }
        }
    }
    /*if (sz == 0)
        return std::numeric_limits<double>::infinity();
    res /= static_cast<double>(sz);*/
    return res;
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror,
            bool isInterp, bool isSubpixel, double subScale)
{
    static const double inf = std::numeric_limits<double>::infinity();
	if (isSubpixel)
	{
		srcImage = resize(srcImage, subScale, isInterp);
	}
    Image out(srcImage.n_rows/3, srcImage.n_cols);
    uint r, g, b;
    double minres1 = inf, minres2 = inf, res;
    int i1_res = 0, j1_res = 0, i2_res = 0, j2_res = 0;
    std::vector<Image> v;
    while (std::min(srcImage.n_rows, srcImage.n_cols) > 400)
    {
    	v.push_back(srcImage);
    	srcImage = resize(srcImage, 0.5, isInterp);
    }
    Image B = srcImage.submatrix(0, 0, srcImage.n_rows/3, srcImage.n_cols);
    Image G = srcImage.submatrix(srcImage.n_rows/3, 0, srcImage.n_rows/3, srcImage.n_cols);
    Image R = srcImage.submatrix(srcImage.n_rows/3 * 2, 0, srcImage.n_rows/3, srcImage.n_cols);
    /*for (int i = -15; i < 15; i++)
    {
        for (int j = -15; j < 15; j++)
        {
            res = MSE(R, G, i, j);
            if (res < minres1)
            {
                minres1 = res; i1_res = i; j1_res = j;
            }
            res = MSE(B, G, i, j);
            if (res < minres2)
            {
                minres2 = res; i2_res = i; j2_res = j;
            }

        }
    }*/
    minres1 = MSE(R, G, 0, 0);
    minres2 = MSE(B, G, 0, 0);
    
    

    bool  keepR = true, keepB = true;
    
    for (int rad = 1; rad < 16; rad++)
    {
    	if (!(keepR || keepB))
    		break;
    	bool changedR = false, changedB = false;
    	for (int l = 0; (l < 2 * rad + 1) && (keepB || keepR); l++)
    	{
    		if (keepR)
    		{
    			res = MSE(R, G, 0 - rad + l, 0 + rad);
    			if (res < minres1)
    			{
    				changedR = true; minres1 = res; i1_res = 0 - rad + l; j1_res = 0 + rad;
    			}
    			res = MSE(R, G, 0 - rad + l, 0 - rad);
    			if (res < minres1)
    			{
    				changedR = true; minres1 = res; i1_res = 0 - rad + l; j1_res = 0 - rad;
    			}
    			res = MSE(R, G, 0 - rad, 0 - rad + l);
    			if (res < minres1)
    			{
    				changedR = true; minres1 = res; i1_res = 0 - rad; j1_res = 0 - rad + l;
    			}
    			res = MSE(R, G, 0 + rad, 0 - rad + l);
    			if (res < minres1)
    			{
    				changedR = true; minres1 = res; i1_res = 0 + rad; j1_res = 0 - rad + l;
    			}
    		}
    		if (keepB)
    		{
    			res = MSE(B, G, 0 - rad + l, 0 + rad);
    			if (res < minres2)
    			{
    				changedB = true; minres2 = res; i2_res = 0 - rad + l; j2_res = 0 + rad;
    			}
    			res = MSE(B, G, 0 - rad + l, 0 - rad);
    			if (res < minres2)
    			{
    				changedB = true; minres2 = res; i2_res = 0 - rad + l; j2_res = 0 - rad;
    			}
    			res = MSE(B, G, 0 - rad, 0 - rad + l);
    			if (res < minres2)
    			{
    				changedB = true; minres2 = res; i2_res = 0 - rad; j2_res = 0 - rad + l;
    			}
    			res = MSE(B, G, 0 + rad, 0 - rad + l);
    			if (res < minres2)
    			{
    				changedB = true; minres2 = res; i2_res = 0 + rad; j2_res = 0 - rad + l;
    			}
    		}
    	}
    	if (!changedR) keepR = false;
    	if (!changedB) keepB = false;
    }
    while (!v.empty())
    {

    	srcImage = v[v.size() - 1];
    	v.pop_back();
    	B = srcImage.submatrix(0, 0, srcImage.n_rows/3, srcImage.n_cols);
	    G = srcImage.submatrix(srcImage.n_rows/3, 0, srcImage.n_rows/3, srcImage.n_cols);
	    R = srcImage.submatrix(srcImage.n_rows/3 * 2, 0, srcImage.n_rows/3, srcImage.n_cols);
	    int i1 = i1_res * 2, i2 = i2_res * 2, j1 = j1_res * 2, j2 = j2_res * 2;
	    i1_res = i1; i2_res = i2; j1_res = j1; j2_res = j2;
        /*minres1 = inf; minres2 = inf;
        for (int i = -2; i < 3; i++)
        {
            for (int j = -2; j < 3; j++)
            {
                res = MSE(R, G, i1 + i, j1 + j);
                if (res < minres1)
                {
                    minres1 = res; i1_res = i1 + i; j1_res = j1 + j;
                }
                res = MSE(B, G, i2 + i, j2 + j);
                if (res < minres2)
                {
                    minres2 = res; i2_res = i2 + i; j2_res = j2 + j;
                }

            }
        }*/
	    minres1 = MSE(R, G, i1, j1);
    	minres2 = MSE(B, G, i2, j2);
    	
	    keepR = true, keepB = true;
	    for (int rad = 1; rad < 2; rad++)
	    {
	    	if (!(keepR || keepB))
	    		break;
	    	bool changedR = false, changedB = false;
	    	for (int l = 0; (l < 2 * rad + 1) && (keepB || keepR); l++)
	    	{
	    		if (keepR)
	    		{
	    			res = MSE(R, G, i1 - rad + l, j1 + rad);
	    			if (res < minres1)
	    			{
	    				changedR = true; minres1 = res; i1_res = i1 - rad + l; j1_res = j1 + rad;
	    			}
	    			res = MSE(R, G, i1 - rad + l, j1 - rad);
	    			if (res < minres1)
	    			{
	    				changedR = true; minres1 = res; i1_res = i1 - rad + l; j1_res = j1 - rad;
	    			}
	    			res = MSE(R, G, i1 - rad, j1 - rad + l);
	    			if (res < minres1)
	    			{
	    				changedR = true; minres1 = res; i1_res = i1 - rad; j1_res = j1 - rad + l;
	    			}
	    			res = MSE(R, G, i1 + rad, j1 - rad + l);
	    			if (res < minres1)
	    			{
	    				changedR = true; minres1 = res; i1_res = i1 + rad; j1_res = j1 - rad + l;
	    			}
	    		}
	    		if (keepB)
	    		{
	    			res = MSE(B, G, i2 - rad + l, j2 + rad);
	    			if (res < minres2)
	    			{
	    				changedB = true; minres2 = res; i2_res = i2 - rad + l; j2_res = j2 + rad;
	    			}
	    			res = MSE(B, G, i2 - rad + l, j2 - rad);
	    			if (res < minres2)
	    			{
	    				changedB = true; minres2 = res; i2_res = i2 - rad + l; j2_res = j2 - rad;
	    			}
	    			res = MSE(B, G, i2 - rad, j2 - rad + l);
	    			if (res < minres2)
	    			{
	    				changedB = true; minres2 = res; i2_res = i2 - rad; j2_res = j2 - rad + l;
	    			}
	    			res = MSE(B, G, i2 + rad, j2 - rad + l);
	    			if (res < minres2)
	    			{
	    				changedB = true; minres2 = res; i2_res = i2 + rad; j2_res = j2 - rad + l;
	    			}
	    		}
	    	}
	    	if (!changedR) keepR = false;
	    	if (!changedB) keepB = false;
	    }
	}
    for (int i = 0; static_cast<unsigned>(i) < srcImage.n_rows/3; i++)
    {
        for (int j = 0; static_cast<unsigned>(j) < srcImage.n_cols; j++)
        {
            r = g = b = 0;
            if ((i + i1_res >= 0) && (static_cast<unsigned>(i + i1_res) < srcImage.n_rows/3) &&
                (j + j1_res >= 0) && (static_cast<unsigned>(j + j1_res) < srcImage.n_cols))
            {
                r = std::get<0>(R(i + i1_res, j + j1_res));
            }
            g = std::get<1>(G(i, j));
            if ((i + i2_res >= 0) && (static_cast<unsigned>(i + i2_res) < srcImage.n_rows/3) &&
                (j + j2_res >= 0) && (static_cast<unsigned>(j + j2_res) < srcImage.n_cols))
            {
                b = std::get<2>(B(i + i2_res, j + j2_res));
            }
            out(i, j) = std::make_tuple(r, g, b);
        }
    }
    if (isSubpixel)
    {
        out = resize(out, 1.0/subScale, isInterp);
    }
    if (isPostprocessing)
    {
        if (postprocessingType == "--gray-world")
            out = gray_world(out);
        else if (postprocessingType == "--unsharp")
        {
            if (isMirror)
            {
                out = mirror(out, 1);
            }
            out = unsharp(out);
            if (isMirror)
            {   
                out = out.submatrix(1, 1, out.n_rows - 2, out.n_cols - 2);
            }
        }
        else if (postprocessingType == "--autocontrast")
        	out = autocontrast(out, fraction);
    }
cout << i1_res << ' ' << j1_res << '\n';
cout << i2_res << ' ' << j2_res << '\n';
    return out;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    Image tmp(src_image.n_rows, src_image.n_cols);
    for (uint i = 1; i < src_image.n_rows - 1; ++i) {
        for (uint j = 1; j < src_image.n_cols - 1; ++j) {
            //auto neighbourhood = submatrix(i - radius, j - radius, size, size);
            double r = 0, g = 0, b = 0;
            for (int k = -1; k < 1 + 1; k++)
            {
                for (int l = - 1; l < 1 + 1; l++)
                {
                    double r1, g1, b1;
                    std::tie(r1, g1, b1) = src_image(i + k, j + l);
                    r += r1*kernel(k + 1, l + 1);
                    g += g1*kernel(k + 1, l + 1);
                    b += b1*kernel(k + 1, l + 1);
                }
            }
            //r /= norm;
            //g /= norm;
            //b /= norm;
            //std::cout << r << ' ' << g << ' ' << b;
            //r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            //g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            //b = (b > 255) ? 255 : (b < 0) ? 0 : b;
            tmp (i,j) = std::make_tuple(r, g, b);
        }
    }
    return tmp;
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    Image tmp(src_image.n_rows, src_image.n_cols);
    for (uint i = 1; i < src_image.n_rows - 1; ++i) {
        for (uint j = 1; j < src_image.n_cols - 1; ++j) {
            //auto neighbourhood = submatrix(i - radius, j - radius, size, size);
            double r = 0, g = 0, b = 0;
            for (int k = -1; k < 1 + 1; k++)
            {
                for (int l = - 1; l < 1 + 1; l++)
                {
                    double r1, g1, b1;
                    std::tie(r1, g1, b1) = src_image(i + k, j + l);
                    r += r1*kernel(k + 1, l + 1);
                    g += g1*kernel(k + 1, l + 1);
                    b += b1*kernel(k + 1, l + 1);
                }
            }
            //r /= norm;
            //g /= norm;
            //b /= norm;
            //std::cout << r << ' ' << g << ' ' << b;
            //r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            //g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            //b = (b > 255) ? 255 : (b < 0) ? 0 : b;
            tmp (i,j) = std::make_tuple(r, g, b);
        }
    }
    return tmp;
}

void sobel_xm(Image src_image, Matrix <int> &m) {
    Matrix<int> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    //Matrix <int> tmp(src_image.n_rows, src_image.n_cols);
    for (uint i = 1; i < src_image.n_rows - 1; ++i) {
        for (uint j = 1; j < src_image.n_cols - 1; ++j) {
            //auto neighbourhood = submatrix(i - radius, j - radius, size, size);
            int r = 0;
            for (int k = -1; k < 1 + 1; k++)
            {
                for (int l = - 1; l < 1 + 1; l++)
                {
                    r += (std::get<0>(src_image(i+k, j+l)))*kernel(k + 1, l + 1);
                }
            }
            m (i,j) = r;
        }
    }
}

void sobel_ym(Image src_image, Matrix <int> &m) {
    Matrix<int> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    //Matrix <int> tmp(src_image.n_rows, src_image.n_cols);
    for (uint i = 1; i < src_image.n_rows - 1; ++i) {
        for (uint j = 1; j < src_image.n_cols - 1; ++j) {
            //auto neighbourhood = submatrix(i - radius, j - radius, size, size);
            int r = 0;
            for (int k = -1; k < 1 + 1; k++)
            {
                for (int l = - 1; l < 1 + 1; l++)
                {
                    r += (std::get<0>(src_image(i+k, j+l)))*kernel(k + 1, l + 1);
                }
            }
            //r /= norm;
            //g /= norm;
            //b /= norm;
            //std::cout << r << ' ' << g << ' ' << b;
            //r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            //g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            //b = (b > 255) ? 255 : (b < 0) ? 0 : b;
            m (i,j) = r;
        }
    }
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {{-1, -4, -1},
                             {-4, 26, -4},
                             {-1, -4, -1}};
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
    double S = 0, Sr = 0, Sg = 0, Sb = 0;                                                                                                                                                                                                                                                                                                                                                                                                      
    for (uint i = 0; i < src_image.n_rows; i++)
    {
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            Sr += std::get<0>(src_image(i,j));
            Sg += std::get<1>(src_image(i,j));
            Sb += std::get<2>(src_image(i,j));
        }
    }
    Sr /= static_cast<double>(src_image.n_rows * src_image.n_cols);
    Sg /= static_cast<double>(src_image.n_rows * src_image.n_cols);
    Sb /= static_cast<double>(src_image.n_rows * src_image.n_cols);
    S = (Sr + Sg + Sb) / 3.0;
    int r,g,b;
    for (uint i = 0; i < src_image.n_rows; i++)
    {
        for (uint j = 0; j < src_image.n_cols; j++)
        {
            std::tie(r, g, b) = src_image(i, j);
            r *= S/Sr;
            g *= S/Sg;
            b *= S/Sb;
            r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            b = (b > 255) ? 255 : (b < 0) ? 0 : b;
            src_image(i,j) = std::make_tuple(r, g, b);
        }
    }
    return src_image;
}

Image resize(Image src_image, double scale, bool isBicubic) {
    double alpha_rows = (scale * src_image.n_rows - 1) / (src_image.n_rows - 1);
    double alpha_cols = (scale * src_image.n_cols - 1) / (src_image.n_cols - 1);
    Image tmp(static_cast<int>(scale * src_image.n_rows), static_cast<int>(scale * src_image.n_cols));
    if (tmp.n_rows * tmp.n_cols == 0)
    {
        cout << "Image became too small!\n";
        return src_image;
    }
    if(!isBicubic)
    {
        int i1, i2, j1, j2, r, g, b;
        double r1, g1, b1;
        src_image = mirror(src_image, 1);
        for (uint i = 0; i < tmp.n_rows; i++)
        {
            for (uint j = 0; j < tmp.n_cols; j++)
            {
                double i0 = static_cast<double>(i) / alpha_rows;
                double j0 = static_cast<double>(j) / alpha_cols;
                i1 = static_cast<int>(i0); //i2 = std::ceil(i0);
                j1 = static_cast<int>(j0); //j2 = std::ceil(j0);
                //if (i1 == i2) i2++;
                //if (j1 == j2) j2++;
                i2 = i1 + 1; j2 = j1 + 1;
                /*if (i == j)
                {
                    cout << i << ' ' << j << ":\n";
                    cout << i0 << ' ' << i1 << ' ' << i2 << '\n';
                    cout << j0 << ' ' << j1 << ' ' << j2 << '\n';
                }*/
                std::tie(r, g, b) = src_image(i1 + 1, j1 + 1);
                r1  = (i2 - i0) * (j2 - j0) * r;
                g1  = (i2 - i0) * (j2 - j0) * g;
                b1  = (i2 - i0) * (j2 - j0) * b;
                std::tie(r, g, b) = src_image(i2 + 1, j1 + 1);
                r1 += (i0 - i1) * (j2 - j0) * r;
                g1 += (i0 - i1) * (j2 - j0) * g;
                b1 += (i0 - i1) * (j2 - j0) * b;
                std::tie(r, g, b) = src_image(i1 + 1, j2 + 1);
                r1 += (i2 - i0) * (j0 - j1) * r;
                g1 += (i2 - i0) * (j0 - j1) * g;
                b1 += (i2 - i0) * (j0 - j1) * b;
                std::tie(r, g, b) = src_image(i2 + 1, j2 + 1);
                r1 += (i0 - i1) * (j0 - j1) * r;
                g1 += (i0 - i1) * (j0 - j1) * g;
                b1 += (i0 - i1) * (j0 - j1) * b;
                tmp(i, j) = std::make_tuple(static_cast <int> (r1), static_cast <int> (g1), static_cast <int> (b1));
             
            }
        }
    }
    else
    {
        int i1, i2, i3, i4, j1, j2, j3, j4, r, g, b;
        double r1, g1, b1;
        src_image = mirror(src_image, 2);
        for (uint i = 0; i < tmp.n_rows; i++)
        {
            for (uint j = 0; j < tmp.n_cols; j++)
            {
                double i0 = static_cast<double>(i) / alpha_rows;
                double j0 = static_cast<double>(j) / alpha_cols;
                //cout << i << '=' << i0 << ' ' << j << '=' << j0 << '\n';
                i2 = static_cast<int>(i0); //i3 = std::ceil(i0);
                j2 = static_cast<int>(j0); //j3 = std::ceil(j0);
                i3 = i2 + 1; j3 = j2 + 1;
                i1 = i2 - 1; j1 = j2 - 1;
                i4 = i3 + 1; j4 = j3 + 1;
                //https://ru.wikipedia.org/wiki/%D0%91%D0%B8%D0%BA%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B0%D1%8F_%D0%B8%D0%BD%D1%82%D0%B5%D1%80%D0%BF%D0%BE%D0%BB%D1%8F%D1%86%D0%B8%D1%8F
                std::tie(r, g, b) = src_image(i2 + 2, j2 + 2);
                r1  =  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*r)/((i2-i1)*(i2-i3)*(i2-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                g1  =  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*g)/((i2-i1)*(i2-i3)*(i2-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                b1  =  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*b)/((i2-i1)*(i2-i3)*(i2-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                std::tie(r, g, b) = src_image(i3 + 2, j2 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*r)/((i3-i1)*(i3-i2)*(i3-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*g)/((i3-i1)*(i3-i2)*(i3-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*b)/((i3-i1)*(i3-i2)*(i3-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                std::tie(r, g, b) = src_image(i2 + 2, j3 + 2);
                r1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*r)/((i2-i1)*(i2-i3)*(i2-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                g1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*g)/((i2-i1)*(i2-i3)*(i2-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                b1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*b)/((i2-i1)*(i2-i3)*(i2-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                std::tie(r, g, b) = src_image(i3 + 2, j3 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*r)/((i3-i1)*(i3-i2)*(i3-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*g)/((i3-i1)*(i3-i2)*(i3-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*b)/((i3-i1)*(i3-i2)*(i3-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                std::tie(r, g, b) = src_image(i1 + 2, j2 + 2);
                r1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*r)/((i1-i2)*(i1-i3)*(i1-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                g1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*g)/((i1-i2)*(i1-i3)*(i1-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                b1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j3)*(j0-j4)*b)/((i1-i2)*(i1-i3)*(i1-i4)*(j2-j1)*(j2-j3)*(j2-j4));
                std::tie(r, g, b) = src_image(i2 + 2, j1 + 2);
                r1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*r)/((i2-i1)*(i2-i3)*(i2-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                g1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*g)/((i2-i1)*(i2-i3)*(i2-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                b1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*b)/((i2-i1)*(i2-i3)*(i2-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                std::tie(r, g, b) = src_image(i1 + 2, j3 + 2);
                r1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*r)/((i1-i2)*(i1-i3)*(i1-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                g1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*g)/((i1-i2)*(i1-i3)*(i1-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                b1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j4)*b)/((i1-i2)*(i1-i3)*(i1-i4)*(j3-j1)*(j3-j2)*(j3-j4));
                std::tie(r, g, b) = src_image(i3 + 2, j1 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*r)/((i3-i1)*(i3-i2)*(i3-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*g)/((i3-i1)*(i3-i2)*(i3-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*b)/((i3-i1)*(i3-i2)*(i3-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                std::tie(r, g, b) = src_image(i2 + 2, j4 + 2);
                r1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*r)/((i2-i1)*(i2-i3)*(i2-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                g1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*g)/((i2-i1)*(i2-i3)*(i2-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                b1 +=  ((i0-i1)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*b)/((i2-i1)*(i2-i3)*(i2-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                std::tie(r, g, b) = src_image(i4 + 2, j2 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j3)*(j0-j4)*r)/((i4-i1)*(i4-i2)*(i4-i3)*(j2-j1)*(j2-j3)*(j2-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j3)*(j0-j4)*g)/((i4-i1)*(i4-i2)*(i4-i3)*(j2-j1)*(j2-j3)*(j2-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j3)*(j0-j4)*b)/((i4-i1)*(i4-i2)*(i4-i3)*(j2-j1)*(j2-j3)*(j2-j4));
                std::tie(r, g, b) = src_image(i1 + 2, j1 + 2);
                r1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*r)/((i1-i2)*(i1-i3)*(i1-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                g1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*g)/((i1-i2)*(i1-i3)*(i1-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                b1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j2)*(j0-j3)*(j0-j4)*b)/((i1-i2)*(i1-i3)*(i1-i4)*(j1-j2)*(j1-j3)*(j1-j4));
                std::tie(r, g, b) = src_image(i4 + 2, j3 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j4)*r)/((i4-i1)*(i4-i2)*(i4-i3)*(j3-j1)*(j3-j2)*(j3-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j4)*g)/((i4-i1)*(i4-i2)*(i4-i3)*(j3-j1)*(j3-j2)*(j3-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j4)*b)/((i4-i1)*(i4-i2)*(i4-i3)*(j3-j1)*(j3-j2)*(j3-j4));
                std::tie(r, g, b) = src_image(i3 + 2, j4 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*r)/((i3-i1)*(i3-i2)*(i3-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*g)/((i3-i1)*(i3-i2)*(i3-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*b)/((i3-i1)*(i3-i2)*(i3-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                std::tie(r, g, b) = src_image(i4 + 2, j1 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j2)*(j0-j3)*(j0-j4)*r)/((i4-i1)*(i4-i2)*(i4-i3)*(j1-j2)*(j1-j3)*(j1-j4));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j2)*(j0-j3)*(j0-j4)*g)/((i4-i1)*(i4-i2)*(i4-i3)*(j1-j2)*(j1-j3)*(j1-j4));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j2)*(j0-j3)*(j0-j4)*b)/((i4-i1)*(i4-i2)*(i4-i3)*(j1-j2)*(j1-j3)*(j1-j4));
                std::tie(r, g, b) = src_image(i1 + 2, j4 + 2);
                r1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*r)/((i1-i2)*(i1-i3)*(i1-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                g1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*g)/((i1-i2)*(i1-i3)*(i1-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                b1 +=  ((i0-i2)*(i0-i3)*(i0-i4)*(j0-j1)*(j0-j2)*(j0-j3)*b)/((i1-i2)*(i1-i3)*(i1-i4)*(j4-j1)*(j4-j2)*(j4-j3));
                std::tie(r, g, b) = src_image(i4 + 2, j4 + 2);
                r1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j3)*r)/((i4-i1)*(i4-i2)*(i4-i3)*(j4-j1)*(j4-j2)*(j4-j3));
                g1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j3)*g)/((i4-i1)*(i4-i2)*(i4-i3)*(j4-j1)*(j4-j2)*(j4-j3));
                b1 +=  ((i0-i1)*(i0-i2)*(i0-i3)*(j0-j1)*(j0-j2)*(j0-j3)*b)/((i4-i1)*(i4-i2)*(i4-i3)*(j4-j1)*(j4-j2)*(j4-j3));
                tmp(i, j) = std::make_tuple(static_cast <int> (r1), static_cast <int> (g1), static_cast <int> (b1));        
            }
        }
    }
    return tmp;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.

    Image tmp(src_image.n_rows, src_image.n_cols);

    const int radiusx = kernel.n_rows / 2;
    const int radiusy = kernel.n_cols / 2;
    const auto sizex = 2 * radiusx + 1;
    const auto sizey = 2 * radiusy + 1;

    const auto start_i = radiusx;
    const auto end_i = src_image.n_rows - radiusx;
    const auto start_j = radiusy;
    const auto end_j = src_image.n_cols - radiusy;
    double norm = 0;
    //std::cout << "NOPE\n" << kernel;
    for (int i = 0; i < sizex; i++)
    {
        for (int j = 0; j < sizey; j++)
        {
            norm += kernel(i, j);
        }
    }
    //std::cout << "CUSTOM NORM: " << norm << '\n';

    for (uint i = start_i; i < end_i; ++i) {
        for (uint j = start_j; j < end_j; ++j) {
            //auto neighbourhood = submatrix(i - radius, j - radius, size, size);
            double r = 0, g = 0, b = 0;
            for (int k = -radiusx; k < radiusx + 1; k++)
            {
                for (int l = - radiusy; l < radiusy + 1; l++)
                {
                    double r1, g1, b1;
                    std::tie(r1, g1, b1) = src_image(i + k, j + l);
                    r += r1*kernel(k + radiusx, l + radiusy);
                    g += g1*kernel(k + radiusx, l + radiusy);
                    b += b1*kernel(k + radiusx, l + radiusy);
                }
            }
            r /= norm;
            g /= norm;
            b /= norm;
            //std::cout << r << ' ' << g << ' ' << b;
            r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            b = (b > 255) ? 255 : (b < 0) ? 0 : b;
            tmp (i,j) = std::make_tuple(r, g, b);
        }
    }
    return tmp;
}

Image autocontrast(Image src_image, double fraction) {
	std::vector<int> v(256);
	int r, g, b, y;
	for (uint i = 0; i < src_image.n_rows; i++)
	{
		for (uint j = 0; j < src_image.n_cols; j++)
		{
			std::tie(r, g, b) = src_image(i, j);
			y = 0.2125*r + 0.7154*g + 0.0721*b;
			//cout << i << ' ' << j << ' ' << y << '\n';
			v[y]++;
		}
	}
	int amount = fraction * (src_image.n_rows*src_image.n_cols) + 1;
	//std::cout << amount << '\n';
	int c1 = 0, c2 = 0, ymin = 255, ymax = 0, minf = 1, maxf = 1;
	for (uint i = 0; i < 256; i++)
	{
		c1 += v[i];
		c2 += v[255 - i];
		if ((c1 >= amount) && minf)
		{
			ymin = (c1 == amount) ? i + 1 : i;
			minf = 0;
		}
		if ((c2 >= amount) && maxf)
		{
			ymax = (c2 == amount) ? 255 - i - 1 : 255 - i;
			maxf = 0;
		}
		if (!(maxf || minf))
			break;
	}
	if (ymax <= ymin)
		return src_image;
	std::cout << ymin << ' ' << ymax << '\n';
	for (uint i = 0; i < src_image.n_rows; i++)
	{
		for (uint j = 0; j < src_image.n_cols; j++)
		{
			std::tie(r, g, b) = src_image(i, j);
			r = ((r - ymin) * 255)/(ymax - ymin);
			g = ((g - ymin) * 255)/(ymax - ymin);
			b = ((b - ymin) * 255)/(ymax - ymin);
			r = (r > 255) ? 255 : (r < 0) ? 0 : r;
            g = (g > 255) ? 255 : (g < 0) ? 0 : g;
            b = (b > 255) ? 255 : (b < 0) ? 0 : b;
			src_image(i, j) = std::make_tuple(r, g, b);
		}
	} 
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    Matrix <double> kernel(2 * radius + 1, 2 * radius + 1);
    for (int i = 0; i < radius + 1; i++)
    {
    	for (int j = 0; j < radius + 1; j++)
    	{
    		double koef = std::exp(-(i*i + j*j)/(2*sigma*sigma));
    		kernel(radius - i, radius - j) = koef;
    		kernel(radius - i, radius + j) = koef;
    		kernel(radius + i, radius - j) = koef;
    		kernel(radius + i, radius + j) = koef;
    	}
    }
    return custom(src_image, kernel);
}
/*
const int maxRange = 15;
int icentral = 0, jcentral = 0, threshold = 15;
*/
double toY(const std::tuple<int, int, int> &t)
{
	return 0.2125*std::get<0>(t) + 0.7154*std::get<1>(t) + 0.0721*std::get<2>(t);
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
	Matrix <double> kernel1(2 * radius + 1, 1), kernel2(1, 2 * radius + 1);
	//double norm = 0;
	for (int i = 0; i < radius + 1; i++)
	{
		double koef = std::exp(-(i*i)/(2*sigma*sigma));
		kernel1(radius - i, 0) = koef;	
		kernel1(radius + i, 0) = koef;
		kernel2(0, radius - i) = koef;
		kernel2(0, radius + i) = koef;
	//	norm += (i == 0) ? koef : 2 * koef;
	}
	src_image = custom(src_image, kernel1);
    return custom(src_image, kernel2);
}

bool med_comp(const std::tuple<int, int, int> &arg1, const std::tuple<int, int, int> &arg2)
{
	return toY(arg1) < toY(arg2);
}

int med(const std::vector <int> &v, const int size)
{
    int sum = 0, end = size * size / 2 + 1;
    for (int i = 0; i < 256; i++)
    {
        sum += v[i];
        if (sum >= end)
            return i;
    }
    return  0;
}

Image median(Image src_image, int radius) {
    Image tmp(src_image.n_rows, src_image.n_cols);

    const int size = 2 * radius + 1;
    
    const int start_i = radius;
    const int end_i = src_image.n_rows - radius;
    const int start_j = radius;
    const int end_j = src_image.n_cols - radius;

    for (int i = start_i; i < end_i; i++)
    {
    	for (int j = start_j; j < end_j; j++)
    	{
            //std::vector <std::tuple <int,int,int>> v;
            std::vector <int> v1(256), v2(256), v3(256);
            int r, g, b;
    		for (int k = 0; k < size; k++)
    		{
    			for (int l = 0; l < size; l++)
    			{
                    std::tie(r, g, b) = src_image(i - radius + k, j - radius + l);
                    v1[r]++; v2[g]++; v3[b]++;
    				//v.push_back(src_image(i - radius + k, j - radius + l));
    				//std::cout << v[radius] << '\n';
    			}
    		}
            //std::sort(v1.begin(), v1.end());
            //std::sort(v2.begin(), v2.end());
            //std::sort(v3.begin(), v3.end());
            tmp(i, j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
			//std::sort(v.begin(), v.end(), med_comp);
			//std::cout << v << ' ' << v2 << ' ' << v3 << '\n';
			//tmp(i, j) = v[size * size / 2];
    	}
    }
    return tmp;
}



Image median_linear(Image src_image, int radius) {
    Image tmp(src_image.n_rows, src_image.n_cols);

    const int size = 2 * radius + 1;
    
    const int start_i = radius;
    const int end_i = src_image.n_rows - radius;
    const int start_j = radius;
    const int end_j = src_image.n_cols - radius;
    
    for (int i = start_i; i < end_i; i++)
    {
        std::vector <int> v1(256), v2(256), v3(256);
        int r, g, b;
        for (int k = i - radius; k <= i + radius; k++)
        {
            for (int l = 0; l < 2 * radius + 1; l++)
            {
                std::tie(r, g, b) = src_image(k, l);
                v1[r]++; v2[g]++; v3[b]++;
            }
        }
        tmp(i, start_j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
        for (int j = start_j + 1; j < end_j; j++)
        {
            for (int k = i - radius; k <= i + radius; k++)
            {
                std::tie(r, g, b) = src_image(k, j - radius - 1);
                v1[r]--; v2[g]--; v3[b]--;
                std::tie(r, g, b) = src_image(k, j + radius);
                v1[r]++; v2[g]++; v3[b]++;
            }
            tmp(i, j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
        }
    }
    return tmp;
}

Image median_const(Image src_image, int radius) {
    Image tmp(src_image.n_rows, src_image.n_cols);

    const int size = 2 * radius + 1;
    
    const int start_i = radius;
    const int end_i = src_image.n_rows - radius;
    const int start_j = radius;
    const int end_j = src_image.n_cols - radius;
    
    //std::vector<std::vector<int>> v1f, v2f, v3f;
    std::vector<std::vector<int>> v1f, v2f, v3f;
    v1f.resize(src_image.n_cols);
    v2f.resize(src_image.n_cols);
    v3f.resize(src_image.n_cols);
    int r, g, b;
{
    std::vector<int> v1(256), v2(256), v3(256);
    for (int j = 0; j < size; j++)
    {
        v1f[j].resize(256);
        v2f[j].resize(256);
        v3f[j].resize(256);
        for (int i = 0; i < size; i++)
        {
            std::tie(r, g, b) = src_image(i, j);
            v1f[j][r]++; v2f[j][g]++; v3f[j][b]++;
        }
        for (int l = 0; l < 256; l++)
        {
            v1[l] += v1f[j][l];
            v2[l] += v2f[j][l];
            v3[l] += v3f[j][l];
        }
    }
    tmp(start_i, start_j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));

    //std::vector<int> v1(256), v2(256), v3(256);
    for (int j = start_j + 1; j < end_j; j++)
    {
        v1f[j + radius].resize(256);
        v2f[j + radius].resize(256);
        v3f[j + radius].resize(256);
        for (int i = 0; i < size; i++)
        {
            std::tie(r, g, b) = src_image(i, j + radius);
            v1f[j + radius][r]++; v2f[j + radius][g]++; v3f[j + radius][b]++;
        }
        for (int l = 0; l < 256; l++)
        {
            v1[l] -= v1f[j - radius - 1][l];
            v2[l] -= v2f[j - radius - 1][l];
            v3[l] -= v3f[j - radius - 1][l];
            v1[l] += v1f[j + radius][l];
            v2[l] += v2f[j + radius][l];
            v3[l] += v3f[j + radius][l];
        }
        tmp(start_i, j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
    }
}
    for (int i = start_i + 1; i < end_i; i++)
    {
        std::vector<int> v1(256), v2(256), v3(256);
        for (int j = 0; j < size; j++)
        {
            std::tie(r, g, b) = src_image(i - radius - 1, j);
            v1f[j][r]--; v2f[j][g]--; v3f[j][b]--;
            std::tie(r, g, b) = src_image(i + radius, j);
            v1f[j][r]++; v2f[j][g]++; v3f[j][b]++;
            for (int l = 0; l < 256; l++)
            {
                v1[l] += v1f[j][l];
                v2[l] += v2f[j][l];
                v3[l] += v3f[j][l];
            }
        }
        tmp(i, start_j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
        for (int j = start_j + 1; j < end_j; j++)
        {
            std::tie(r, g, b) = src_image(i - radius - 1, j + radius);
            v1f[j + radius][r]--; v2f[j + radius][g]--; v3f[j + radius][b]--;
            std::tie(r, g, b) = src_image(i + radius, j + radius);
            v1f[j + radius][r]++; v2f[j + radius][g]++; v3f[j + radius][b]++;       
            for (int l = 0; l < 256; l++)
            {
                v1[l] -= v1f[j - radius - 1][l];
                v2[l] -= v2f[j - radius - 1][l];
                v3[l] -= v3f[j - radius - 1][l];
                v1[l] += v1f[j + radius][l];
                v2[l] += v2f[j + radius][l];
                v3[l] += v3f[j + radius][l];
            }   
            tmp(i, j) = std::make_tuple(med(v1, size), med(v2, size), med(v3, size));
        }
    }

    return tmp;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    Image tmpi = mirror(src_image, 4);
	//save_image(src_image, "S1.bmp");   
    tmpi = gaussian_separable(tmpi, 1.4, 2);
    //save_image(tmpi, "S2.bmp");
    //cout << "@@\n";
    Matrix<int> X(tmpi.n_rows, tmpi.n_cols), Y(tmpi.n_rows, tmpi.n_cols);
    sobel_xm(tmpi, X); sobel_ym(tmpi, Y);
    //save_image(X, "S3.bmp"); save_image(Y, "S4.bmp");
    Matrix <double> tmp(tmpi.n_rows, tmpi.n_cols), tan(tmpi.n_rows, tmpi.n_cols), G(tmpi.n_rows, tmpi.n_cols);
    for (uint i = 0; i < tmpi.n_rows; i++)
    {
    	for (uint j = 0; j < tmpi.n_cols; j++)
    	{
    		double Ix = (X(i, j)), Iy = (Y(i, j));
    		tmp(i, j) = std::sqrt(Ix * Ix + Iy * Iy);
    		tan(i, j) = std::atan2(Iy, Ix) / M_PI;
    	}
    }
    //cout << tmp;
    for (uint i = 2; i < tmpi.n_rows - 2; i++)
    {
    	for (uint j = 2; j < tmpi.n_cols - 2; j++)
    	{
            double t = tan(i, j), g = tmp(i, j);
    		if(((t >= -1.0) && (t < -7.0/8)) || ((t >= 7.0/8) && (t <= 1.0)) || ((t >= -1.0/8) && (t < 1.0/8)))
    		{
                if ((g > tmp(i, j + 1)) && (g > tmp (i, j - 1))) G(i, j) = 1;
                {
                    if (g > threshold2) G(i, j) = 1;
                    else if (g > threshold1) G(i, j) = 0.5;
                    else G(i, j) = 0;
                }
    			//G(i, j) = ((g > tmp(i, j + 1)) && (g > tmp (i, j - 1))) ? (g > threshold2) ? 1 : (g > threshold1) ? 0.5 : 0 : 0;
    		}
    		else if (((t >= 1.0/8) && (t < -3.0/8)) || ((t >= -7.0/8) && (t < -5.0/8)))
    		{
                if ((g > tmp(i, j + 1)) && (g > tmp (i, j - 1))) G(i, j) = 1;
                {
                    if (g > threshold2) G(i, j) = 1;
                    else if (g > threshold1) G(i, j) = 0.5;
                    else G(i, j) = 0;
                }
    			//G(i, j) = ((g > tmp(i + 1, j - 1)) && (g > tmp(i - 1, j + 1))) ? (g > threshold2) ? 1 : (g > threshold1) ? 0.5 : 0 : 0;
    		}
    		else if (((t >= -5.0/8) && (t < -3.0/8)) || ((t >= 3.0/8) && (t < 5.0/8)))
    		{
                if ((g > tmp(i - 1, j)) && (g > tmp(i + 1, j))) G(i, j) = 1;
                {
                    if (g > threshold2) G(i, j) = 1;
                    else if (g > threshold1) G(i, j) = 0.5;
                    else G(i, j) = 0;
                }
    			//G(i, j) = ((g > tmp(i - 1, j)) && (g > tmp(i + 1, j))) ? (g > threshold2) ? 1 : (g > threshold1) ? 0.5 : 0 : 0;
    		}
    		else if (((t >= -3.0/8) && (t < -1.0/8)) || ((t >= 5.0/8) && (t < 7.0/8)))
    		{
                if ((g > tmp(i + 1, j + 1)) && (g > tmp(i - 1, j - 1))) G(i, j) = 1;
                {
                    if (g > threshold2) G(i, j) = 1;
                    else if (g > threshold1) G(i, j) = 0.5;
                    else G(i, j) = 0;
                }
    			//G(i, j) = ((g > tmp(i + 1, j + 1)) && (g > tmp(i - 1, j - 1))) ? (g > threshold2) ? 1 : (g > threshold1) ? 0.5 : 0 : 0;
    		}
    	}
    }
    /*Image out(tmpi.n_rows - 8, tmpi.n_cols - 8);
    for (uint i = 0; i < out.n_rows; i++)
    {
    	for (uint j = 0; j < out.n_cols; j++)
    	{
            //out(i, j) = std::make_tuple((G(i + 1, j + 1) > 0.75) ? 255 : 0, ((G(i + 1, j + 1) > 0.25) && (G(i + 1,j + 1) < 0.75)) ? 255 : 0, 0);
    		out(i, j) = std::make_tuple(((G(i + 4, j + 4) > 0.25) && (G(i + 4, j + 4) < 0.75)) ? 255 : 0, (G(i + 4, j + 4) > 0.75) ? 255 : 0, 0);
    	}
    }*/
      //  cout << "1231\n";
    G = G.submatrix(4, 4, G.n_rows - 8, G.n_cols - 8);
    std::vector<std::pair<uint, uint>> vl(G.n_rows), vu(G.n_cols), vr(G.n_rows), vd(G.n_cols);
    //int counter;
    for (uint i = 0; i < G.n_rows; i++)
    {
        for ( vl[i].first = 0;  (vl[i].first < G.n_cols / 20) && (G(i, vl[i].first) < 0.25);  vl[i].first++);
        for (;  (vl[i].first < G.n_cols / 20) && (G(i, vl[i].first) >= 0.75); vl[i].first++);
        for ( vl[i].second = vl[i].first;  (vl[i].second < G.n_cols / 20) && (G(i, vl[i].second) < 0.25); vl[i].second++);
        for (;  (vl[i].second < G.n_cols / 20) && (G(i, vl[i].second) >= 0.75); vl[i].second++);
        for ( vr[i].first = 0;  (vr[i].first <  G.n_cols / 20) && (G(i, G.n_cols - 1 - vr[i].first) < 0.25); vr[i].first++);
        for (;  (vr[i].first < G.n_cols / 20) && (G(i, G.n_cols - 1 - vr[i].first) >= 0.75); vr[i].first++);
        for ( vr[i].second = vr[i].first;  (vr[i].second < G.n_cols / 20) && (G(i, G.n_cols - 1 - vr[i].second) < 0.25); vr[i].second++);
        for (;  (vr[i].second < G.n_cols / 20) && (G(i, G.n_cols - 1 - vr[i].second) >= 0.75); vr[i].second++);
    }
    for (uint i = 0; i < G.n_cols; i++)
    {
        for ( vu[i].first = 0;  (vu[i].first < G.n_rows / 20) && (G(vu[i].first, i) < 0.25);  vu[i].first++);
        for (;  (vu[i].first < G.n_rows / 20) && (G(vu[i].first, i) >= 0.75); vu[i].first++);
        for ( vu[i].second = vu[i].first;  (vu[i].second < G.n_rows / 20) && (G(vu[i].second, i) < 0.25); vu[i].second++);
        for (;  (vu[i].second < G.n_rows / 20) && (G(vu[i].second, i) >= 0.75); vu[i].second++);
        for ( vd[i].first = 0;  (vd[i].first <  G.n_rows / 20) && (G(G.n_rows - 1 - vd[i].first, i) < 0.25); vd[i].first++);
        for (;  (vd[i].first < G.n_rows / 20) && (G(G.n_rows - 1 - vd[i].first, i) >= 0.75); vd[i].first++);
        for ( vd[i].second = vd[i].first;  (vd[i].second < G.n_rows / 20) && (G(G.n_rows - 1 - vd[i].second, i) < 0.25); vd[i].second++);
        for (;  (vd[i].second < G.n_rows / 20) && (G(G.n_rows - 1 - vd[i].second, i) >= 0.75); vd[i].second++);
    }
    std::sort(vl.begin(), vl.end());
    std::sort(vr.begin(), vr.end());
    std::sort(vu.begin(), vu.end());
    std::sort(vd.begin(), vd.end());
    //cout << "1260\n";
    Image out = src_image.submatrix(4 + vu[vu.size()/2].second, 4 + vl[vl.size()/2].second, G.n_rows - vu[vu.size()/2].second - vd[vd.size()/2].second, G.n_cols - vl[vl.size()/2].second - vr[vr.size()/2].second);

    return out;
}
