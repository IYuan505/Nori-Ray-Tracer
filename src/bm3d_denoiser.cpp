#include <nori/color.h>
#include <nori/vector.h>
#include <nori/denoiser.h>

NORI_NAMESPACE_BEGIN

std::vector<float> kaiserWindow = {
    0.1924f, 0.2989f, 0.3846f, 0.4325f, 0.4325f, 0.3845f, 0.2989f, 0.1924f,
    0.2989f, 0.4642f, 0.5974f, 0.6717f, 0.6717f, 0.5974f, 0.4642f, 0.2989f,
    0.3846f, 0.5974f, 0.7688f, 0.8644f, 0.8644f, 0.7689f, 0.5974f, 0.3846f,
    0.4325f, 0.6717f, 0.8644f, 0.9718f, 0.9718f, 0.8644f, 0.6717f, 0.4325f,
    0.4325f, 0.6717f, 0.8644f, 0.9718f, 0.9718f, 0.8644f, 0.6717f, 0.4325f,
    0.3846f, 0.5974f, 0.7688f, 0.8644f, 0.8644f, 0.7689f, 0.5974f, 0.3846f,
    0.2989f, 0.4642f, 0.5974f, 0.6717f, 0.6717f, 0.5974f, 0.4642f, 0.2989f,
    0.1924f, 0.2989f, 0.3846f, 0.4325f, 0.4325f, 0.3845f, 0.2989f, 0.1924f
};

bool ComparaisonFirst(std::pair<float,unsigned> pair1, std::pair<float,unsigned> pair2) {
    return pair1.first < pair2.first;
}

int closest_power_of_2(unsigned n) {
    unsigned r = 1;
    while (r * 2 <= n)
        r *= 2;
    return r;
}

unsigned log2(unsigned N) {
    unsigned k = 1;
    unsigned n = 0;
    while (k < N) {
        k *= 2;
        n++;
    }
    return n;
}

void color_space_transform(std::vector<float> &img, unsigned width, unsigned height, unsigned chnls, bool rgb2yuv) {
    std::vector<float> tmp;
    tmp.resize(chnls * width * height);
    unsigned red   = 0;
    unsigned green = width * height;
    unsigned blue  = width * height * 2;
    if (rgb2yuv) {
        for (unsigned k = 0; k < width * height; k++){
            tmp[k + red  ] =  0.299f   * img[k + red] + 0.587f   * img[k + green] + 0.114f   * img[k + blue];
            tmp[k + green] = -0.14713f * img[k + red] - 0.28886f * img[k + green] + 0.436f   * img[k + blue];
            tmp[k + blue ] =  0.615f   * img[k + red] - 0.51498f * img[k + green] - 0.10001f * img[k + blue];
        }
    } else {
        for (unsigned k = 0; k < width * height; k++) {
            tmp[k + red  ] = img[k + red] + 1.13983f * img[k + blue];
            tmp[k + green] = img[k + red] - 0.39465f * img[k + green] - 0.5806f * img[k + blue];
            tmp[k + blue ] = img[k + red] + 2.03211f * img[k + green];
        }
    }
    for (unsigned k = 0; k < width * height * chnls; k++)
        img[k] = tmp[k];
}

/**
 * \brief Denoise the rendered image
 *
 * The denoise class provides denoising function for the bitmap
 */
class BM3D_denoiser: public Denoiser {
public:
    BM3D_denoiser(const PropertyList &propList) {}

    void load_image(Bitmap *noisy, std::vector<float> &img, unsigned *width, unsigned *height, unsigned *chnls) const {
        *width = noisy->cols();
        *height = noisy->rows();
        *chnls = 3;
        unsigned tot_size = (*width) * (*height) * (*chnls);
        img.resize(tot_size);
        unsigned red   = 0;
        unsigned green = (*width) * (*height);
        unsigned blue  = (*width) * (*height) * 2;
        for (unsigned i = 0; i < *height; ++i) {
            for (unsigned j = 0; j < *width; ++j) {
                Color3f tonemapped = noisy->coeffRef(i, j).toSRGB();
                img[i * (*width) + j + red] = clamp(255.f * tonemapped[0], 0.f, 255.f);
                img[i * (*width) + j + green] = clamp(255.f * tonemapped[1], 0.f, 255.f);
                img[i * (*width) + j + blue] = clamp(255.f * tonemapped[2], 0.f, 255.f);
            }
        }
    }

    Bitmap *store_bitmap(std::vector<float> img, unsigned width, unsigned height) const {
        Bitmap *results = new Bitmap(Vector2i(width, height));
        unsigned red   = 0;
        unsigned green = width * height;
        unsigned blue  = width * height * 2;
        for (unsigned i = 0; i < height; ++i) {
            for (unsigned j = 0; j < width; ++j) {
                results->coeffRef(i, j)[0] = img[i * width + j + red];
                results->coeffRef(i, j)[1] = img[i * width + j + green];
                results->coeffRef(i, j)[2] = img[i * width + j + blue];
            }
        }
        return results;
    }

    void symetrize(std::vector<float> &img, std::vector<float> &img_sym, unsigned width, unsigned height, 
        unsigned chnls, unsigned N) const {
        unsigned w = width + 2 * N;
        unsigned h = height + 2 * N;

        img_sym.resize(w * h * chnls);

        for (unsigned c = 0; c < chnls; c++)
        {
            unsigned dc = c * width * height;
            unsigned dc_2 = c * w * h + N * w + N;

            // Center of the image
            for (unsigned i = 0; i < height; i++)
                for (unsigned j = 0; j < width; j++, dc++)
                    img_sym[dc_2 + i * w + j] = img[dc];

            // Top and bottom
            dc_2 = c * w * h;
            for (unsigned j = 0; j < w; j++, dc_2++)
                for (unsigned i = 0; i < N; i++)
                {
                    img_sym[dc_2 + i * w] = img_sym[dc_2 + (2 * N - i - 1) * w];
                    img_sym[dc_2 + (h - i - 1) * w] = img_sym[dc_2 + (h - 2 * N + i) * w];
                }

            // Right and left
            dc_2 = c * w * h;
            for (unsigned i = 0; i < h; i++)
            {
                unsigned di = dc_2 + i * w;
                for (unsigned j = 0; j < N; j++)
                {
                    img_sym[di + j] = img_sym[di + 2 * N - j - 1];
                    img_sym[di + w - j - 1] = img_sym[di + w - 2 * N + j];
                }
            }
        }
    }

    void estimate_sigma(float sigma, std::vector<float> &sigma_table, unsigned chnls) const {
        // map for YUV color space
        sigma_table[0] = sqrtf(0.299f * 0.299f + 0.587f * 0.587f + 0.114f * 0.114f) * sigma;
        sigma_table[1] = sqrtf(0.14713f * 0.14713f + 0.28886f * 0.28886f + 0.436f * 0.436f) * sigma;
        sigma_table[2] = sqrtf(0.615f * 0.615f + 0.51498f * 0.51498f + 0.10001f * 0.10001f) * sigma;
    }

    void ind_initialize(std::vector<unsigned> &ind_set, unsigned max_size, unsigned N, unsigned step) const  {
        ind_set.clear();
        unsigned ind = N;
        while (ind < max_size - N)
        {
            ind_set.push_back(ind);
            ind += step;
        }
        if (ind_set.back() < max_size - N - 1)
            ind_set.push_back(max_size - N - 1);
    }

    void bior15_coef(std::vector<float> &lp1, std::vector<float> &hp1, std::vector<float> &lp2, std::vector<float> &hp2) const {
        float coef_norm = 1.f / (sqrtf(2.f) * 128.f);
        float sqrt2_inv = 1.f / sqrtf(2.f);

        lp1.resize(10); 
        lp1[0] = 3.f;   lp1[1] = -3.f; lp1[2] = -22.f; lp1[3] = 22.f; lp1[4] = 128.f;
        lp1[5] = 128.f; lp1[6] = 22.f; lp1[7] = -22.f; lp1[8] = -3.f; lp1[9] = 3.f;

        hp1.resize(10);
        hp1[0] = 0.f;       hp1[1] = 0.f; hp1[2] = 0.f; hp1[3] = 0.f; hp1[4] = -sqrt2_inv;
        hp1[5] = sqrt2_inv; hp1[6] = 0.f; hp1[7] = 0.f; hp1[8] = 0.f; hp1[9] = 0.f;

        lp2.resize(10);     
        lp2[0] = 0.f;       lp2[1] = 0.f; lp2[2] = 0.f; lp2[3] = 0.f; lp2[4] = sqrt2_inv;
        lp2[5] = sqrt2_inv; lp2[6] = 0.f; lp2[7] = 0.f; lp2[8] = 0.f; lp2[9] = 0.f;

        hp2.resize(10);
        hp2[0] = 3.f;    hp2[1] = 3.f;  hp2[2] = -22.f; hp2[3] = -22.f; hp2[4] = 128.f;
        hp2[5] = -128.f; hp2[6] = 22.f; hp2[7] = 22.f;  hp2[8] = -3.f;  hp2[9] = -3.f;

        for (unsigned k = 0; k < 10; k++)
        {
            lp1[k] *= coef_norm;
            hp2[k] *= coef_norm;
        }
    }

    void precompute_BM(std::vector<std::vector<unsigned> > &patch_table, std::vector<float> &img, unsigned width, unsigned height,
        unsigned kHW, unsigned NHW, unsigned nHW, unsigned pHW, float tauMatch) const {
        //! Declarations
        unsigned Ns = 2 * nHW + 1;
        float threshold = tauMatch * kHW * kHW;
        std::vector<float> diff_table(width * height);
        std::vector<std::vector<float> > sum_table((nHW + 1) * Ns, std::vector<float> (width * height, 2 * threshold));
        if (patch_table.size() != width * height)
            patch_table.resize(width * height);
        std::vector<unsigned> row_ind;
        ind_initialize(row_ind, height - kHW + 1, nHW, pHW);
        std::vector<unsigned> column_ind;
        ind_initialize(column_ind, width - kHW + 1, nHW, pHW);

        //! For each possible distance, precompute inter-patches distance
        for (unsigned di = 0; di <= nHW; di++) {
            for (unsigned dj = 0; dj < Ns; dj++) {
                int dk = (int) (di * width + dj) - (int) nHW;
                unsigned ddk = di * Ns + dj;

                //! Process the image containing the square distance between pixels
                for (unsigned i = nHW; i < height - nHW; i++) {
                    unsigned k = i * width + nHW;
                    for (unsigned j = nHW; j < width - nHW; j++, k++)
                        diff_table[k] = (img[k + dk] - img[k]) * (img[k + dk] - img[k]);
                }

                //! Compute the sum for each patches, using the method of the integral images
                unsigned dn = nHW * width + nHW;
                //! 1st patch, top left corner
                float value = 0.0f;
                for (unsigned p = 0; p < kHW; p++) {
                    unsigned pq = p * width + dn;
                    for (unsigned q = 0; q < kHW; q++, pq++)
                        value += diff_table[pq];
                }
                sum_table[ddk][dn] = value;

                //! 1st row, top
                for (unsigned j = nHW + 1; j < width - nHW; j++) {
                    unsigned ind = nHW * width + j - 1;
                    float sum = sum_table[ddk][ind];
                    for (unsigned p = 0; p < kHW; p++)
                        sum += diff_table[ind + p * width + kHW] - diff_table[ind + p * width];
                    sum_table[ddk][ind + 1] = sum;
                }

                //! General case
                for (unsigned i = nHW + 1; i < height - nHW; i++) {
                    unsigned ind = (i - 1) * width + nHW;
                    float sum = sum_table[ddk][ind];
                    //! 1st column, left
                    for (unsigned q = 0; q < kHW; q++)
                        sum += diff_table[ind + kHW * width + q] - diff_table[ind + q];
                    sum_table[ddk][ind + width] = sum;

                    //! Other columns
                    unsigned k = i * width + nHW + 1;
                    unsigned pq = (i + kHW - 1) * width + kHW - 1 + nHW + 1;
                    for (unsigned j = nHW + 1; j < width - nHW; j++, k++, pq++) {
                        sum_table[ddk][k] =
                              sum_table[ddk][k - 1]
                            + sum_table[ddk][k - width]
                            - sum_table[ddk][k - 1 - width]
                            + diff_table[pq]
                            - diff_table[pq - kHW]
                            - diff_table[pq - kHW * width]
                            + diff_table[pq - kHW - kHW * width];
                    }

                }
            }
        }

        //! Precompute Bloc Matching
        std::vector<std::pair<float, unsigned> > table_distance;
        //! To avoid reallocation
        table_distance.reserve(Ns * Ns);

        for(unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++) {
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
                //! Initialization
                unsigned k_r = row_ind[ind_i] * width + column_ind[ind_j];
                table_distance.clear();
                patch_table[k_r].clear();

                //! Threshold distances in order to keep similar patches
                for (int dj = -(int) nHW; dj <= (int) nHW; dj++)
                {
                    for (int di = 0; di <= (int) nHW; di++)
                        if (sum_table[dj + nHW + di * Ns][k_r] < threshold)
                            table_distance.push_back(std::make_pair(
                                        sum_table[dj + nHW + di * Ns][k_r]
                                      , k_r + di * width + dj));

                    for (int di = - (int) nHW; di < 0; di++)
                        if (sum_table[-dj + nHW + (-di) * Ns][k_r + di * width + dj] < threshold)
                            table_distance.push_back(std::make_pair(
                                        sum_table[-dj + nHW + (-di) * Ns][k_r + di * width + dj]
                                      , k_r + di * width + dj));
                }

                // We need a power of 2 for the number of similar patches,
                // because of the Welsh-Hadamard transform on the third dimension.
                // We assume that NHW is already a power of 2
                unsigned nSx_r = (NHW > table_distance.size() ?
                                        closest_power_of_2(table_distance.size()) : NHW);

                // To avoid problem
                if (nSx_r == 1 && table_distance.size() == 0) {
                    table_distance.push_back(std::make_pair(0, k_r));
                }

                // Sort patches according to their distance to the reference one
                partial_sort(table_distance.begin(), table_distance.begin() + nSx_r,
                                                table_distance.end(), ComparaisonFirst);

                // Keep a maximum of NHW similar patches
                for (unsigned n = 0; n < nSx_r; n++)
                    patch_table[k_r].push_back(table_distance[n].second);

            }
        }
    }

    void per_ext_ind(std::vector<unsigned> &ind_per, unsigned N, unsigned L) const{
        for (unsigned k = 0; k < N; k++)
            ind_per[k + L] = k;

        int ind1 = (N - L);
        while (ind1 < 0)
            ind1 += N;
        unsigned ind2 = 0;
        unsigned k = 0;
        while(k < L) {
            ind_per[k] = (unsigned) ind1;
            ind_per[k + L + N] = ind2;
            ind1 = ((unsigned) ind1 < N - 1 ? (unsigned) ind1 + 1 : 0);
            ind2 = (ind2 < N - 1 ? ind2 + 1 : 0);
            k++;
        }
    }

    void bior_2d_forward(std::vector<float> &input, std::vector<float> &output, unsigned N, unsigned d_i,
        unsigned r_i, unsigned d_o, std::vector<float> &lpd, std::vector<float> &hpd) const {
        // Initializing output
        for (unsigned i = 0; i < N; i++)
            for (unsigned j = 0; j < N; j++)
                output[i * N + j + d_o] = input[i * r_i + j + d_i];

        unsigned iter_max = log2(N);
        unsigned N_1 = N;
        unsigned N_2 = N / 2;
        unsigned S_1 = lpd.size();
        unsigned S_2 = S_1 / 2 - 1;

        for (unsigned iter = 0; iter < iter_max; iter++)
        {
            // Periodic extension index initialization
            std::vector<float> tmp(N_1 + 2 * S_2);
            std::vector<unsigned> ind_per(N_1 + 2 * S_2);
            per_ext_ind(ind_per, N_1, S_2);

            // Implementing row filtering
            for (unsigned i = 0; i < N_1; i++) {
                // Periodic extension of the signal in row
                for (unsigned j = 0; j < tmp.size(); j++)
                    tmp[j] = output[d_o + i * N + ind_per[j]];

                // Low and High frequencies filtering
                for (unsigned j = 0; j < N_2; j++) {
                    float v_l = 0.0f, v_h = 0.0f;
                    for (unsigned k = 0; k < S_1; k++)
                    {
                        v_l += tmp[k + j * 2] * lpd[k];
                        v_h += tmp[k + j * 2] * hpd[k];
                    }
                    output[d_o + i * N + j] = v_l;
                    output[d_o + i * N + j + N_2] = v_h;
                }
            }

            // Implementing column filtering
            for (unsigned j = 0; j < N_1; j++) {
                // Periodic extension of the signal in column
                for (unsigned i = 0; i < tmp.size(); i++)
                    tmp[i] = output[d_o + j + ind_per[i] * N];

                // Low and High frequencies filtering
                for (unsigned i = 0; i < N_2; i++) {
                    float v_l = 0.0f, v_h = 0.0f;
                    for (unsigned k = 0; k < S_1; k++) {
                        v_l += tmp[k + i * 2] * lpd[k];
                        v_h += tmp[k + i * 2] * hpd[k];
                    }
                    output[d_o + j + i * N] = v_l;
                    output[d_o + j + (i + N_2) * N] = v_h;
                }
            }

            // Sizes update
            N_1 /= 2;
            N_2 /= 2;
        }
    }

    void bior_2d_inverse(std::vector<float> &signal, unsigned N, unsigned d_s, 
        std::vector<float> &lpr, std::vector<float> &hpr) const {
        unsigned iter_max = log2(N);
        unsigned N_1 = 2;
        unsigned N_2 = 1;
        unsigned S_1 = lpr.size();
        unsigned S_2 = S_1 / 2 - 1;

        for (unsigned iter = 0; iter < iter_max; iter++) {
            std::vector<float> tmp(N_1 + S_2 * N_1);
            std::vector<unsigned> ind_per(N_1 + 2 * S_2 * N_2);
            per_ext_ind(ind_per, N_1, S_2 * N_2);

            // Implementing column filtering
            for (unsigned j = 0; j < N_1; j++) {
                // Periodic extension of the signal in column
                for (unsigned i = 0; i < tmp.size(); i++)
                    tmp[i] = signal[d_s + j + ind_per[i] * N];

                // Low and High frequencies filtering
                for (unsigned i = 0; i < N_2; i++) {
                    float v_l = 0.0f, v_h = 0.0f;
                    for (unsigned k = 0; k < S_1; k++) {
                        v_l += lpr[k] * tmp[k * N_2 + i];
                        v_h += hpr[k] * tmp[k * N_2 + i];
                    }

                    signal[d_s + i * 2 * N + j] = v_h;
                    signal[d_s + (i * 2 + 1) * N + j] = v_l;
                }
            }

            // Implementing row filtering
            for (unsigned i = 0; i < N_1; i++) {
                // Periodic extension of the signal in row
                for (unsigned j = 0; j < tmp.size(); j++)
                    tmp[j] = signal[d_s + i * N + ind_per[j]];

                // Low and High frequencies filtering
                for (unsigned j = 0; j < N_2; j++) {
                    float v_l = 0.0f, v_h = 0.0f;
                    for (unsigned k = 0; k < S_1; k++) {
                        v_l += lpr[k] * tmp[k * N_2 + j];
                        v_h += hpr[k] * tmp[k * N_2 + j];
                    }

                    signal[d_s + i * N + j * 2] = v_h;
                    signal[d_s + i * N + j * 2 + 1] = v_l;
                }
            }

            // Sizes update
            N_1 *= 2;
            N_2 *= 2;
        }
    }

    void bior_2d_inverse_warpper(std::vector<float> &group_3D_table, unsigned kHW, std::vector<float> &lpr, 
        std::vector<float> &hpr) const {
        unsigned kHW_2 = kHW * kHW;
        unsigned N = group_3D_table.size() / kHW_2;

        // Bior process
        for (unsigned n = 0; n < N; n++)
            bior_2d_inverse(group_3D_table, kHW, n * kHW_2, lpr, hpr);
    }

    void bior_2d_process(std::vector<float> &bior_table_2D, std::vector<float> &img, unsigned nHW, unsigned width,
        unsigned height, unsigned chnls, unsigned kHW, unsigned i_r, unsigned step, unsigned i_min, unsigned i_max,
        std::vector<float> &lpd, std::vector<float> &hpd) const {
        unsigned kHW_2 = kHW * kHW;

        // If i_r == ns, then we have to process all Bior1.5 transforms
        if (i_r == i_min || i_r == i_max) {
            for (unsigned c = 0; c < chnls; c++)
            {
                unsigned dc = c * width * height;
                unsigned dc_p = c * kHW_2 * width * (2 * nHW + 1);
                for (unsigned i = 0; i < 2 * nHW + 1; i++)
                    for (unsigned j = 0; j < width - kHW; j++)
                    {
                        bior_2d_forward(img, bior_table_2D, kHW, dc +
                                 (i_r + i - nHW) * width + j, width,
                                 dc_p + (i * width + j) * kHW_2, lpd, hpd);
                    }
            }
        } else {
            unsigned ds = step * width * kHW_2;

            // Re-use of Bior1.5 already processed
            for (unsigned c = 0; c < chnls; c++) {
                unsigned dc = c * width * (2 * nHW + 1) * kHW_2;
                for (unsigned i = 0; i < 2 * nHW + 1 - step; i++)
                    for (unsigned j = 0; j < width - kHW; j++)
                        for (unsigned k = 0; k < kHW_2; k++)
                            bior_table_2D[k + (i * width + j) * kHW_2 + dc] =
                                bior_table_2D[k + (i * width + j) * kHW_2 + dc + ds];
            }

            // Compute the new Bior
            for (unsigned c = 0; c < chnls; c++) {
                unsigned dc   = c * width * height;
                unsigned dc_p = c * kHW_2 * width * (2 * nHW + 1);
                for (unsigned i = 0; i < step; i++)
                    for (unsigned j = 0; j < width - kHW; j++) {
                        bior_2d_forward(img, bior_table_2D, kHW, dc +
                                 (i + 2 * nHW + 1 - step + i_r - nHW) * width + j,
                                 width, dc_p + ((i + 2 * nHW + 1 - step)
                                 * width + j) * kHW_2, lpd, hpd);
                    }
            }
        }
    }

    void hadamard_transform(std::vector<float> &vec, std::vector<float> &tmp, unsigned N, unsigned D) const {
        if (N == 1)
            return;
        else if (N == 2) {
            float a = vec[D + 0];
            float b = vec[D + 1];
            vec[D + 0] = a + b;
            vec[D + 1] = a - b;
        } else {
            unsigned n = N / 2;
            for (unsigned k = 0; k < n; k++) {
                float a = vec[D + 2 * k];
                float b = vec[D + 2 * k + 1];
                vec[D + k] = a + b;
                tmp[k] = a - b;
            }
            for (unsigned k = 0; k < n; k++)
                vec[D + n + k] = tmp[k];

            hadamard_transform(vec, tmp, n, D);
            hadamard_transform(vec, tmp, n, D + n);
        }
    }

    void ht_filtering_hadamard(std::vector<float> &group_3D, std::vector<float> &tmp, unsigned nSx_r, unsigned kHard,
        unsigned chnls, std::vector<float> const& sigma_table, float lambdaHard3D, std::vector<float> &weight_table) const {
        unsigned kHard_2 = kHard * kHard;
        for (unsigned c = 0; c < chnls; c++)
            weight_table[c] = 0.0f;
        float coef_norm = sqrtf((float) nSx_r);
        float coef = 1.0f / (float) nSx_r;

        // Process the Welsh-Hadamard transform on the 3rd dimension
        for (unsigned n = 0; n < kHard_2 * chnls; n++)
            hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

        // Hard Thresholding
        for (unsigned c = 0; c < chnls; c++) {
            unsigned dc = c * nSx_r * kHard_2;
            float T = lambdaHard3D * sigma_table[c] * coef_norm;
            for (unsigned k = 0; k < kHard_2 * nSx_r; k++) {
                if (k < 1 || fabs(group_3D[k + dc]) > T)
                    weight_table[c]++;
                else
                    group_3D[k + dc] = 0.0f;
            }
        }

        //! Process of the Welsh-Hadamard inverse transform
        for (unsigned n = 0; n < kHard_2 * chnls; n++)
            hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

        for (unsigned k = 0; k < group_3D.size(); k++)
            group_3D[k] *= coef;
    }

    void wiener_filtering_hadamard(std::vector<float> &group_3D_img, std::vector<float> &group_3D_est, std::vector<float> &tmp,
        unsigned nSx_r, unsigned kWien, unsigned chnls, std::vector<float> &sigma_table, std::vector<float> &weight_table) const {
        unsigned kWien_2 = kWien * kWien;
        float coef = 1.0f / (float) nSx_r;

        for (unsigned c = 0; c < chnls; c++)
            weight_table[c] = 0.0f;

        // Process the Welsh-Hadamard transform on the 3rd dimension
        for (unsigned n = 0; n < kWien_2 * chnls; n++)
        {
            hadamard_transform(group_3D_img, tmp, nSx_r, n * nSx_r);
            hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);
        }

        // Wiener Filtering
        for (unsigned c = 0; c < chnls; c++)
        {
            unsigned dc = c * nSx_r * kWien_2;
            group_3D_est[dc] = group_3D_img[dc] * coef;
            // Add the weight corresponding to the DC components that were not passed through the Wiener filter
            weight_table[c] += 1; 
            for (unsigned k = 1; k < kWien_2 * nSx_r; k++) {
                float value = group_3D_est[dc + k] * group_3D_est[dc + k] * coef;
                value /= (value + sigma_table[c] * sigma_table[c]);
                group_3D_est[k + dc] = group_3D_img[k + dc] * value * coef;
                weight_table[c] += (value*value);
            }
        }

        // Process of the Welsh-Hadamard inverse transform
        for (unsigned n = 0; n < kWien_2 * chnls; n++)
            hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);
    }

    void sd_weighting(std::vector<float> const& group_3D, unsigned nSx_r, unsigned kHW, unsigned chnls,
        std::vector<float> &weight_table) const {
        unsigned N = nSx_r * kHW * kHW;

        for (unsigned c = 0; c < chnls; c++) {
            float mean = 0.0f;
            float std  = 0.0f;

            // Compute the sum and the square sum
            for (unsigned k = 0; k < N; k++) {
                mean += group_3D[k];
                std  += group_3D[k] * group_3D[k];
            }

            // Sample standard deviation (Bessel's correction)
            float res = (std - mean * mean / (float) N) / (float) (N - 1);

            // Return the weight as used in the aggregation
            weight_table[c] = (res > 0.0f ? 1.0f / sqrtf(res) : 0.0f);
        }
    }

    void bm3d_1st_step(float sigma, std::vector<float> & img_noisy, std::vector<float> &img_basic, unsigned width, 
        unsigned height, unsigned chnls, unsigned nHard, unsigned kHard, unsigned NHard, unsigned pHard) const {
        // Estimatation of sigma on each channel
        std::vector<float> sigma_table(chnls);
        estimate_sigma(sigma, sigma_table, chnls);

        // Parameters initialization
        float lambdaHard3D = 2.7f;            // Threshold for Hard Thresholding
        float tauMatch = sigma_table[0] < 35.0f ? 2500 : 5000; // Threshold used to determinate similarity between patches

        // Initialization for convenience
        std::vector<unsigned> row_ind;
        ind_initialize(row_ind, height - kHard + 1, nHard, pHard);
        std::vector<unsigned> column_ind;
        ind_initialize(column_ind, width - kHard + 1, nHard, pHard);
        unsigned kHard_2 = kHard * kHard;
        std::vector<float> group_3D_table(chnls * kHard_2 * NHard * column_ind.size());
        std::vector<float> wx_r_table;
        wx_r_table.reserve(chnls * column_ind.size());
        std::vector<float> hadamard_tmp(NHard);

        // Resize basic
        img_basic.resize(img_noisy.size());

        // Preprocessing of Bior table
        std::vector<float> lpd, hpd, lpr, hpr;
        bior15_coef(lpd, hpd, lpr, hpr);

        // For aggregation part
        std::vector<float> denominator(width * height * chnls, 0.0f);
        std::vector<float> numerator  (width * height * chnls, 0.0f);

        // Precompute Bloc-Matching
        std::vector<std::vector<unsigned> > patch_table;
        precompute_BM(patch_table, img_noisy, width, height, kHard, NHard, nHard, pHard, tauMatch);

        // table_2D[p * N + q + (i * width + j) * kHard_2 + c * (2 * nHard + 1) * width * kHard_2]
        std::vector<float> table_2D((2 * nHard + 1) * width * chnls * kHard_2, 0.0f);

        // Loop on i_r
        for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
        {
            unsigned i_r = row_ind[ind_i];

            bior_2d_process(table_2D, img_noisy, nHard, width, height, chnls,
                kHard, i_r, pHard, row_ind[0], row_ind.back(), lpd, hpd);

            wx_r_table.clear();
            group_3D_table.clear();

            // Loop on j_r
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
            {
                unsigned j_r = column_ind[ind_j];
                unsigned k_r = i_r * width + j_r;

                // Number of similar patches
                unsigned nSx_r = patch_table[k_r].size();

                // Build of the 3D group
                std::vector<float> group_3D(chnls * nSx_r * kHard_2, 0.0f);
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++)
                    {
                        unsigned ind = patch_table[k_r][n] + (nHard - i_r) * width;
                        for (unsigned k = 0; k < kHard_2; k++)
                            group_3D[n + k * nSx_r + c * kHard_2 * nSx_r] =
                                table_2D[k + ind * kHard_2 + c * kHard_2 * (2 * nHard + 1) * width];
                    }

                // HT filtering of the 3D group
                std::vector<float> weight_table(chnls);
                ht_filtering_hadamard(group_3D, hadamard_tmp, nSx_r, kHard, chnls, sigma_table,
                    lambdaHard3D, weight_table);

                sd_weighting(group_3D, nSx_r, kHard, chnls, weight_table);

                // Save the 3D group. The BIOR 2D inverse will be done after.
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++)
                        for (unsigned k = 0; k < kHard_2; k++)
                            group_3D_table.push_back(group_3D[n + k * nSx_r + c * kHard_2 * nSx_r]);

                // Save weighting
                for (unsigned c = 0; c < chnls; c++)
                    wx_r_table.push_back(weight_table[c]);

            } // End of loop on j_r
            
            bior_2d_inverse_warpper(group_3D_table, kHard, lpr, hpr);

            //! Registration of the weighted estimation
            unsigned dec = 0;
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
                unsigned j_r   = column_ind[ind_j];
                unsigned k_r   = i_r * width + j_r;
                unsigned nSx_r = patch_table[k_r].size();
                for (unsigned c = 0; c < chnls; c++) {
                    for (unsigned n = 0; n < nSx_r; n++) {
                        unsigned k = patch_table[k_r][n] + c * width * height;
                        for (unsigned p = 0; p < kHard; p++)
                            for (unsigned q = 0; q < kHard; q++) {
                                unsigned ind = k + p * width + q;
                                numerator[ind] += kaiserWindow[p * kHard + q]
                                                * wx_r_table[c + ind_j * chnls]
                                                * group_3D_table[p * kHard + q + n * kHard_2
                                                      + c * kHard_2 * nSx_r + dec];
                                denominator[ind] += kaiserWindow[p * kHard + q]
                                                * wx_r_table[c + ind_j * chnls];
                            }
                    }
                }
                dec += nSx_r * chnls * kHard_2;
            }

        } // End of loop on i_r

        // Final reconstruction
        for (unsigned k = 0; k < width * height * chnls; k++)
            img_basic[k] = numerator[k] / denominator[k];
    }

    void bm3d_2nd_step(float sigma, std::vector<float> &img_noisy, std::vector<float> &img_basic, std::vector<float> &img_denoised,
        unsigned width, unsigned height, unsigned chnls, unsigned nWien, unsigned kWien, unsigned NWien, unsigned pWien) const {
        // Estimatation of sigma on each channel
        std::vector<float> sigma_table(chnls);
        estimate_sigma(sigma, sigma_table, chnls);

        // Parameters initialization
        float tauMatch = sigma_table[0] < 35.0f ? 400 : 3500; // threshold used to determinate similarity between patches

        // Initialization for convenience
        std::vector<unsigned> row_ind;
        ind_initialize(row_ind, height - kWien + 1, nWien, pWien);
        std::vector<unsigned> column_ind;
        ind_initialize(column_ind, width - kWien + 1, nWien, pWien);
        unsigned kWien_2 = kWien * kWien;
        std::vector<float> group_3D_table(chnls * kWien_2 * NWien * column_ind.size());
        std::vector<float> wx_r_table;
        wx_r_table.reserve(chnls * column_ind.size());
        std::vector<float> tmp(NWien);

        // Resize denoised
        img_denoised.resize(img_noisy.size());

        // Preprocessing of Bior table
        std::vector<float> lpd, hpd, lpr, hpr;
        bior15_coef(lpd, hpd, lpr, hpr);

        // For aggregation part
        std::vector<float> denominator(width * height * chnls, 0.0f);
        std::vector<float> numerator  (width * height * chnls, 0.0f);

        // Precompute Bloc-Matching
        std::vector<std::vector<unsigned> > patch_table;
        precompute_BM(patch_table, img_basic, width, height, kWien, NWien, nWien, pWien, tauMatch);


        // table_2D[p * N + q + (i * width + j) * kWien_2 + c * (2 * ns + 1) * width * kWien_2]
        std::vector<float> table_2D_img((2 * nWien + 1) * width * chnls * kWien_2, 0.0f);
        std::vector<float> table_2D_est((2 * nWien + 1) * width * chnls * kWien_2, 0.0f);

        // Loop on i_r
        for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++) {
            unsigned i_r = row_ind[ind_i];
            
            bior_2d_process(table_2D_img, img_noisy, nWien, width, height,
                chnls, kWien, i_r, pWien, row_ind[0], row_ind.back(), lpd, hpd);
            bior_2d_process(table_2D_est, img_basic, nWien, width, height,
                chnls, kWien, i_r, pWien, row_ind[0], row_ind.back(), lpd, hpd);

            wx_r_table.clear();
            group_3D_table.clear();

            // Loop on j_r
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
                unsigned j_r = column_ind[ind_j];
                unsigned k_r = i_r * width + j_r;

                // Number of similar patches
                unsigned nSx_r = patch_table[k_r].size();

                // Build of the 3D group
                std::vector<float> group_3D_est(chnls * nSx_r * kWien_2, 0.0f);
                std::vector<float> group_3D_img(chnls * nSx_r * kWien_2, 0.0f);
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++) {
                        unsigned ind = patch_table[k_r][n] + (nWien - i_r) * width;
                        for (unsigned k = 0; k < kWien_2; k++)
                        {
                            group_3D_est[n + k * nSx_r + c * kWien_2 * nSx_r] =
                                table_2D_est[k + ind * kWien_2 + c * kWien_2 * (2 * nWien + 1) * width];
                            group_3D_img[n + k * nSx_r + c * kWien_2 * nSx_r] =
                                table_2D_img[k + ind * kWien_2 + c * kWien_2 * (2 * nWien + 1) * width];
                        }
                    }

                // Wiener filtering of the 3D group
                std::vector<float> weight_table(chnls);
                wiener_filtering_hadamard(group_3D_img, group_3D_est, tmp, nSx_r, kWien,
                                                chnls, sigma_table, weight_table);

                // 3D weighting using Standard Deviation
                sd_weighting(group_3D_est, nSx_r, kWien, chnls, weight_table);

                // Save the 3D group. The DCT 2D inverse will be done after.
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++)
                        for (unsigned k = 0; k < kWien_2; k++)
                            group_3D_table.push_back(group_3D_est[n + k * nSx_r + c * kWien_2 * nSx_r]);

                // Save weighting
                for (unsigned c = 0; c < chnls; c++)
                    wx_r_table.push_back(weight_table[c]);

            } // End of loop on j_r

            //  Apply 2D bior inverse
            bior_2d_inverse_warpper(group_3D_table, kWien, lpr, hpr);

            // Registration of the weighted estimation
            unsigned dec = 0;
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
                unsigned j_r   = column_ind[ind_j];
                unsigned k_r   = i_r * width + j_r;
                unsigned nSx_r = patch_table[k_r].size();
                for (unsigned c = 0; c < chnls; c++) {
                    for (unsigned n = 0; n < nSx_r; n++) {
                        unsigned k = patch_table[k_r][n] + c * width * height;
                        for (unsigned p = 0; p < kWien; p++)
                            for (unsigned q = 0; q < kWien; q++)
                            {
                                unsigned ind = k + p * width + q;
                                numerator[ind] += kaiserWindow[p * kWien + q]
                                                * wx_r_table[c + ind_j * chnls]
                                                * group_3D_table[p * kWien + q + n * kWien_2
                                                      + c * kWien_2 * nSx_r + dec];
                                denominator[ind] += kaiserWindow[p * kWien + q]
                                                * wx_r_table[c + ind_j * chnls];
                            }
                    }
                }
                dec += nSx_r * chnls * kWien_2;
            }

        } // End of loop on i_r

        // Final reconstruction
        for (unsigned k = 0; k < width * height * chnls; k++)
            img_denoised[k] = numerator[k] / denominator[k];
    }

    /// Denoise the bitmap
    Bitmap *denoise(Bitmap *noisy, float sigma) const {
        std::vector<float> img_noisy, img_basic, img_denoised;
        unsigned width, height, chnls;
        load_image(noisy, img_noisy, &width, &height, &chnls);
        img_basic.resize(img_noisy.size());
        img_denoised.resize(img_noisy.size());

        // Parameters
        unsigned nHard = 16; //! Half size of the search window
        unsigned nWien = 16; //! Half size of the search window
        unsigned NHard = 16; //! Must be a power of 2
        unsigned NWien = 32; //! Must be a power of 2
        unsigned pHard = 3;
        unsigned pWien = 3;
        unsigned kHard = 8;
        unsigned kWien = 8;

        // Transformation to YUV color space
        color_space_transform(img_noisy, width, height, chnls, true);

        // Add boundaries and symetrize them
        unsigned h_b = height + 2 * nHard;
        unsigned w_b = width  + 2 * nHard;
        std::vector<float> img_sym_noisy, img_sym_basic, img_sym_denoised;
        symetrize(img_noisy, img_sym_noisy, width, height, chnls, nHard);

        // Denoising, 1st Step
        bm3d_1st_step(sigma, img_sym_noisy, img_sym_basic, w_b, h_b, chnls, nHard, kHard, NHard, pHard);

        // To avoid boundaries problem
        for (unsigned c = 0; c < chnls; c++) {
            unsigned dc_b = c * w_b * h_b + nHard * w_b + nHard;
            unsigned dc = c * width * height;
            for (unsigned i = 0; i < height; i++)
                for (unsigned j = 0; j < width; j++, dc++)
                    img_basic[dc] = img_sym_basic[dc_b + i * w_b + j];
        }
        symetrize(img_basic, img_sym_basic, width, height, chnls, nHard);

        bm3d_2nd_step(sigma, img_sym_noisy, img_sym_basic, img_sym_denoised, w_b, h_b, chnls, nWien, kWien, NWien, pWien);

        // Obtention of img_denoised
        for (unsigned c = 0; c < chnls; c++)
        {
            unsigned dc_b = c * w_b * h_b + nWien * w_b + nWien;
            unsigned dc = c * width * height;
            for (unsigned i = 0; i < height; i++)
                for (unsigned j = 0; j < width; j++, dc++)
                    img_denoised[dc] = img_sym_denoised[dc_b + i * w_b + j];
        }

        // Inverse color space transform to RGB
        color_space_transform(img_denoised, width, height, chnls, false);
        color_space_transform(img_noisy, width, height, chnls, false);
        color_space_transform(img_basic, width, height, chnls, false);

        Bitmap *results = store_bitmap(img_denoised, width, height);
        return results;
    }
    std::string toString() const {
        return "BM3D denoiser\n";
    }
};

NORI_REGISTER_CLASS(BM3D_denoiser, "bm3d_denoiser");
NORI_NAMESPACE_END
