#include <nori/color.h>
#include <nori/vector.h>
#include <nori/denoiser.h>

NORI_NAMESPACE_BEGIN

std::vector<float> kaiserWindow = {
    0.6239f, 0.6633f, 0.6981f, 0.7277f, 0.7520f, 0.7704f, 0.7828f, 0.7891f, 0.7891f, 0.7828f, 0.7704f, 0.7520f, 0.7277f, 0.6981f, 0.6633f, 0.6239f, 
    0.6633f, 0.7052f, 0.7422f, 0.7737f, 0.7995f, 0.8191f, 0.8323f, 0.8389f, 0.8389f, 0.8323f, 0.8191f, 0.7995f, 0.7737f, 0.7422f, 0.7052f, 0.6633f, 
    0.6981f, 0.7422f, 0.7811f, 0.8143f, 0.8414f, 0.8620f, 0.8759f, 0.8829f, 0.8829f, 0.8759f, 0.8620f, 0.8414f, 0.8143f, 0.7811f, 0.7422f, 0.6981f, 
    0.7277f, 0.7737f, 0.8143f, 0.8489f, 0.8772f, 0.8987f, 0.9132f, 0.9205f, 0.9205f, 0.9132f, 0.8987f, 0.8772f, 0.8489f, 0.8143f, 0.7737f, 0.7277f, 
    0.7520f, 0.7995f, 0.8414f, 0.8772f, 0.9064f, 0.9286f, 0.9435f, 0.9511f, 0.9511f, 0.9435f, 0.9286f, 0.9064f, 0.8772f, 0.8414f, 0.7995f, 0.7520f, 
    0.7704f, 0.8191f, 0.8620f, 0.8987f, 0.9286f, 0.9513f, 0.9667f, 0.9744f, 0.9744f, 0.9667f, 0.9513f, 0.9286f, 0.8987f, 0.8620f, 0.8191f, 0.7704f, 
    0.7828f, 0.8323f, 0.8759f, 0.9132f, 0.9435f, 0.9667f, 0.9823f, 0.9901f, 0.9901f, 0.9823f, 0.9667f, 0.9435f, 0.9132f, 0.8759f, 0.8323f, 0.7828f, 
    0.7891f, 0.8389f, 0.8829f, 0.9205f, 0.9511f, 0.9744f, 0.9901f, 0.9980f, 0.9980f, 0.9901f, 0.9744f, 0.9511f, 0.9205f, 0.8829f, 0.8389f, 0.7891f, 
    0.7891f, 0.8389f, 0.8829f, 0.9205f, 0.9511f, 0.9744f, 0.9901f, 0.9980f, 0.9980f, 0.9901f, 0.9744f, 0.9511f, 0.9205f, 0.8829f, 0.8389f, 0.7891f, 
    0.7828f, 0.8323f, 0.8759f, 0.9132f, 0.9435f, 0.9667f, 0.9823f, 0.9901f, 0.9901f, 0.9823f, 0.9667f, 0.9435f, 0.9132f, 0.8759f, 0.8323f, 0.7828f, 
    0.7704f, 0.8191f, 0.8620f, 0.8987f, 0.9286f, 0.9513f, 0.9667f, 0.9744f, 0.9744f, 0.9667f, 0.9513f, 0.9286f, 0.8987f, 0.8620f, 0.8191f, 0.7704f, 
    0.7520f, 0.7995f, 0.8414f, 0.8772f, 0.9064f, 0.9286f, 0.9435f, 0.9511f, 0.9511f, 0.9435f, 0.9286f, 0.9064f, 0.8772f, 0.8414f, 0.7995f, 0.7520f, 
    0.7277f, 0.7737f, 0.8143f, 0.8489f, 0.8772f, 0.8987f, 0.9132f, 0.9205f, 0.9205f, 0.9132f, 0.8987f, 0.8772f, 0.8489f, 0.8143f, 0.7737f, 0.7277f, 
    0.6981f, 0.7422f, 0.7811f, 0.8143f, 0.8414f, 0.8620f, 0.8759f, 0.8829f, 0.8829f, 0.8759f, 0.8620f, 0.8414f, 0.8143f, 0.7811f, 0.7422f, 0.6981f, 
    0.6633f, 0.7052f, 0.7422f, 0.7737f, 0.7995f, 0.8191f, 0.8323f, 0.8389f, 0.8389f, 0.8323f, 0.8191f, 0.7995f, 0.7737f, 0.7422f, 0.7052f, 0.6633f, 
    0.6239f, 0.6633f, 0.6981f, 0.7277f, 0.7520f, 0.7704f, 0.7828f, 0.7891f, 0.7891f, 0.7828f, 0.7704f, 0.7520f, 0.7277f, 0.6981f, 0.6633f, 0.6239f
};
std::vector<float> LoD = {
    3, -3, -22, 22, 128,
    128, 22, -22, -3, 3
};
std::vector<float> HiD = {
    0.f, 0.f, 0.f, 0.f, -INV_SQRT_TWO,
    INV_SQRT_TWO, 0.f, 0.f, 0.f, 0.f
};
std::vector<float> LoR = {
    0.f, 0.f, 0.f, 0.f, INV_SQRT_TWO,
    INV_SQRT_TWO, 0.f, 0.f, 0.f, 0.f
};
std::vector<float> HiR = {
    3, 3, -22, -22, 128,
    -128, 22, 22, -3, -3
};

bool ComparaisonFirst(std::pair<float,unsigned> pair1, std::pair<float,unsigned> pair2) {
    return pair1.first < pair2.first;
}

int closestPower2(unsigned n) {
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

void indInitialize(std::vector<unsigned> &ind_set, unsigned max_size, unsigned N, unsigned step) {
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

void perExtInd(std::vector<unsigned> &ind_per, unsigned N, unsigned L) {
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

void colorSpaceTransform(std::vector<float> &img, unsigned width, unsigned height, unsigned chnls, bool rgb2opp) {
    std::vector<float> tmp;
    tmp.resize(chnls * width * height);
    unsigned red   = 0;
    unsigned green = width * height;
    unsigned blue  = width * height * 2;
    if (rgb2opp) {
        for (unsigned k = 0; k < width * height; k++){
            tmp[k + red  ] =  0.333f * img[k + red] + 0.333f * img[k + green] + 0.333f * img[k + blue];
            tmp[k + green] =  0.500f * img[k + red] + 0.000f * img[k + green] - 0.500f * img[k + blue];
            tmp[k + blue ] =  0.250f * img[k + red] - 0.500f * img[k + green] + 0.250f * img[k + blue];
        }
    } else {
        for (unsigned k = 0; k < width * height; k++) {
            tmp[k + red  ] = 1.0f * img[k + red] + 1.0f * img[k + green] + 0.666f * img[k + blue];
            tmp[k + green] = 1.0f * img[k + red] + 0.0f * img[k + green] - 1.333f * img[k + blue];
            tmp[k + blue ] = 1.0f * img[k + red] - 1.0f * img[k + green] + 0.666f * img[k + blue];
        }
    }
    for (unsigned k = 0; k < width * height * chnls; k++)
        img[k] = tmp[k];
}

void mapSigma(float sigma, std::vector<float> &sigma_table) {
    // map for OPP color space
    sigma_table[0] = sqrtf(0.333f * 0.333f + 0.333f * 0.333f + 0.333f * 0.333f) * sigma;
    sigma_table[1] = sqrtf(0.500f * 0.500f + 0.000f * 0.000f + 0.500f * 0.500f) * sigma;
    sigma_table[2] = sqrtf(0.250f * 0.250f + 0.500f * 0.500f + 0.250f * 0.250f) * sigma;
}

void symmetricPadding(std::vector<float> &img, std::vector<float> &img_sym, unsigned width, unsigned height, 
    unsigned chnls, unsigned N) {
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

void bior2dForward(std::vector<float> &input, std::vector<float> &output, unsigned N, 
    unsigned d_i, unsigned r_i, unsigned d_o) {
    // Initializing output
    for (unsigned i = 0; i < N; i++)
        for (unsigned j = 0; j < N; j++)
            output[i * N + j + d_o] = input[i * r_i + j + d_i];

    unsigned iter_max = log2(N);
    unsigned N_1 = N;
    unsigned N_2 = N / 2;
    unsigned S_1 = LoD.size();
    unsigned S_2 = S_1 / 2 - 1;

    for (unsigned iter = 0; iter < iter_max; iter++)
    {
        // Periodic extension index initialization
        std::vector<float> tmp(N_1 + 2 * S_2);
        std::vector<unsigned> ind_per(N_1 + 2 * S_2);
        perExtInd(ind_per, N_1, S_2);

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
                    v_l += tmp[k + j * 2] * LoD[k];
                    v_h += tmp[k + j * 2] * HiD[k];
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
                    v_l += tmp[k + i * 2] * LoD[k];
                    v_h += tmp[k + i * 2] * HiD[k];
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

void bior2dInverse(std::vector<float> &signal, unsigned N, unsigned d_s) {
    unsigned iter_max = log2(N);
    unsigned N_1 = 2;
    unsigned N_2 = 1;
    unsigned S_1 = LoR.size();
    unsigned S_2 = S_1 / 2 - 1;

    for (unsigned iter = 0; iter < iter_max; iter++) {
        std::vector<float> tmp(N_1 + S_2 * N_1);
        std::vector<unsigned> ind_per(N_1 + 2 * S_2 * N_2);
        perExtInd(ind_per, N_1, S_2 * N_2);

        // Implementing column filtering
        for (unsigned j = 0; j < N_1; j++) {
            // Periodic extension of the signal in column
            for (unsigned i = 0; i < tmp.size(); i++)
                tmp[i] = signal[d_s + j + ind_per[i] * N];

            // Low and High frequencies filtering
            for (unsigned i = 0; i < N_2; i++) {
                float v_l = 0.0f, v_h = 0.0f;
                for (unsigned k = 0; k < S_1; k++) {
                    v_l += tmp[k * N_2 + i] * LoR[k];
                    v_h += tmp[k * N_2 + i] * HiR[k];
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
                    v_l += tmp[k * N_2 + j] * LoR[k];
                    v_h += tmp[k * N_2 + j] * HiR[k];
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

void hadamardTransform(std::vector<float> &vec, std::vector<float> &tmp, unsigned N, unsigned D) {
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

        hadamardTransform(vec, tmp, n, D);
        hadamardTransform(vec, tmp, n, D + n);
    }
}

/**
 * \brief Denoise the rendered image
 *
 * The denoise class provides denoising function for the bitmap
 */
class BM3D_denoiser: public Denoiser {
public:
    BM3D_denoiser(const PropertyList &propList) {}

    void loadImage(Bitmap *noisy, std::vector<float> &img, unsigned *width, unsigned *height, unsigned *chnls) const {
        unsigned tot_size;         // total size of image array
        unsigned red, green, blue; // start index of red, green, blue channel
        
        *width = noisy->cols();
        *height = noisy->rows();
        *chnls = 3;
        tot_size = (*width) * (*height) * (*chnls);
        red = 0, green = (*width) * (*height), blue = (*width) * (*height) * 2;
        img.resize(tot_size);
        for (unsigned i = 0; i < *height; ++i) {
            for (unsigned j = 0; j < *width; ++j) {
                Color3f tonemapped = noisy->coeffRef(i, j).toSRGB();
                img[red + i * (*width) + j] = clamp(255.f * tonemapped[0], 0.f, 255.f);
                img[green + i * (*width) + j] = clamp(255.f * tonemapped[1], 0.f, 255.f);
                img[blue + i * (*width) + j] = clamp(255.f * tonemapped[2], 0.f, 255.f);
            }
        }
    }

    Bitmap *storeBitmap(std::vector<float> img, unsigned width, unsigned height) const {
        Bitmap *results = new Bitmap(Vector2i(width, height));
        unsigned red   = 0;
        unsigned green = width * height;
        unsigned blue  = width * height * 2;

        for (unsigned i = 0; i < height; ++i) {
            for (unsigned j = 0; j < width; ++j) {
                results->coeffRef(i, j)[0] = img[red + i * width + j];
                results->coeffRef(i, j)[1] = img[green + i * width + j];
                results->coeffRef(i, j)[2] = img[blue + i * width + j];
            }
        }
        return results;
    }

    void drawSamples(std::vector<unsigned> &samples, Sampler *sampler, 
        unsigned width, unsigned height, unsigned chnls) const {
        unsigned nSamples = 1e5;

        samples.resize(nSamples);
        for (unsigned sample = 0; sample < nSamples; sample++)
        {
            unsigned x = unsigned(sampler->next1D() * width * height * chnls) % width;
            unsigned y = unsigned(sampler->next1D() * width * height * chnls) % height;
            unsigned z = unsigned(sampler->next1D() * width * height * chnls) % chnls;
            samples[sample] = z * width * height + y * width + x;;
        }
        std::sort(samples.begin(), samples.end());
    }

    void getPatchStd(std::vector<unsigned> samples, std::vector<float> &patchStd,
        std::vector<float> img, unsigned width, unsigned height, unsigned chnls) const {
        std::vector<float> img_sym;
        unsigned paddingSize = 8;
        unsigned patchSize = 15;
        unsigned patchLength = patchSize*patchSize;
        symmetricPadding(img, img_sym, width, height, chnls, paddingSize);

        patchStd.resize(samples.size());
        for (unsigned sample = 0; sample < samples.size(); sample++) {
            unsigned idx = samples[sample];
            int sample_z = idx / (width * height);
            int sample_y = (idx - sample_z * width * height) / width;
            int sample_x = idx - sample_z * width * height - sample_y * width;

            float sum = 0;
            for (unsigned y = sample_y; y < sample_y + patchSize; y++) {
                unsigned pos = sample_z * (width + 2 * paddingSize) * (height + 2 * paddingSize)
                     + y * (width + 2 * paddingSize) + sample_x;
                for (unsigned x = 0; x < patchSize; x++)
                    sum += img_sym[pos + x];
            }
            float mean = sum / patchLength;
            float variance = 0.0f;
            for (unsigned y = sample_y; y < sample_y + patchSize; y++) {
                unsigned pos = sample_z * (width + 2 * paddingSize) * (height + 2 * paddingSize)
                     + y * (width + 2 * paddingSize) + sample_x;
                for (unsigned x = 0; x < patchSize; x++)
                    variance += (img_sym[pos + x] - mean) * (img_sym[pos + x] - mean);
            }
            variance = variance / patchLength;
            patchStd[sample] = sqrtf(variance);
        }
    }

    void estimatePatchNoise(std::vector<float>patchStd, int sampleCount, float &sigma) const {
        unsigned fs = 0, ls = sampleCount-1;
        unsigned* assignment = (unsigned*) malloc(sampleCount * sizeof(unsigned));
        float minimum = 1e20f;
        float maximum = -1e20f;
        float mean1, mean2;
        float sum1 = 0.0f, sum2 = 0.0f, cnt1 = 0.0f, cnt2 = 0.0f;
        float var1 = 0.0f, var2 = 0.0f;
        bool change = true;
        unsigned iters_taken = 0;
        unsigned maxiters = 1000;

        for (unsigned idx = fs; idx <= ls; idx++) {
            float std = patchStd[idx];
            if (std < minimum) minimum = std;
            if (std > maximum) maximum = std;
        }
        mean1 = minimum + 0.25 * (maximum - minimum);
        mean2 = maximum - 0.25 * (maximum - minimum);

        for (unsigned idx = fs; idx <= ls; idx++) {
            float std = patchStd[idx];
            if (std == 0.0f) assignment[idx] = -1;
            else if (fabs(std - mean1) < fabs(std - mean2)) {
                sum1 += std;
                cnt1++;
                assignment[idx] = 0;
            } else {
                sum2 += std;
                cnt2++;
                assignment[idx] = 1;
            }
        }

        if(cnt1 > 0.0f) mean1 = sum1 / cnt1;
        if(cnt2 > 0.0f) mean2 = sum2 / cnt2;
        if (cnt1 == 0 && cnt2 == 0) {
            sigma = 0.0f; return;
        }
        
        for (unsigned idx = fs; idx <= ls; idx++) {
            float std = patchStd[idx];
            if (assignment[idx] == 0)
                var1 += (std - mean1) * (std - mean1);
            else if (assignment[idx] == 1)
                var2 += (std - mean2) * (std - mean2);
        }
        var1 /= std::max(1.0f , cnt1);
        var2 /= std::max(1.0f , cnt2);

        while (change && var1 > 0.0f && var2 > 0.0f) {
            change = false;
            iters_taken++;
            sum1 = sum2 = cnt1 = cnt2 = 0.0f;

            for (unsigned idx = fs; idx <= ls; idx++) {
                float std = patchStd[idx];
                if (std == 0.0f) continue;

                float prob1 = 1./(sqrtf(var1 * 2 * M_PI)) * expf(-((std - mean1) * (std - mean1) / (2 * var1)));
                float prob2 = 1./(sqrtf(var2 * 2 * M_PI)) * expf(-((std - mean2) * (std - mean2) / (2 * var2)));

                if (prob2 > prob1) {
                    if (assignment[idx] != 1) {
                        change = true;
                        assignment[idx] = 1;
                    }
                    sum2 += std;
                    cnt2++;
                } else {
                    if (assignment[idx] != 0) {
                        change = true;
                        assignment[idx] = 0;
                    }
                    sum1 += std;
                    cnt1++;
                }
            }
            if(cnt1 > 0.0f) mean1 = sum1 / cnt1;
            if(cnt2 > 0.0f) mean2 = sum2 / cnt2;

            var1 = 0.0f;
            var2 = 0.0f;
            for (unsigned idx = fs; idx <= ls; idx++) {
                float std = patchStd[idx];
                if (assignment[idx] == 0)
                    var1 += (std - mean1) * (std - mean1);
                else if (assignment[idx] == 1)
                    var2 += (std - mean2) * (std - mean2);
            }
            var1 /= std::max(1.0f , cnt1);
            var2 /= std::max(1.0f , cnt2);

            if (iters_taken == maxiters) break;
        }
        free(assignment);
        if (mean1 <= mean2)
            sigma = mean1 + 3 * sqrt(var1);
        else 
            sigma = mean2 + 3 * sqrt(var2);
    }

    void noiseEstimation(std::vector<float> img, unsigned width, unsigned height, unsigned chnls,
        Sampler *sampler, float &sigma) const {
        std::vector<unsigned> samples;
        std::vector<float> patchStd;
        drawSamples(samples, sampler, width, height, chnls);
        getPatchStd(samples, patchStd, img, width, height, chnls);
        estimatePatchNoise(patchStd, samples.size(), sigma);
    }

    void precomputeBM(std::vector<std::vector<unsigned> > &patch_table, std::vector<float> &img, unsigned width, unsigned height,
        unsigned kHW, unsigned NHW, unsigned nHW, unsigned pHW, float tauMatch) const {
        unsigned Ns = 2 * nHW + 1;
        float threshold = tauMatch * kHW * kHW;
        std::vector<float> diff_table(width * height);
        std::vector<std::vector<float> > sum_table((nHW + 1) * Ns, std::vector<float> (width * height, 2 * threshold));
        if (patch_table.size() != width * height)
            patch_table.resize(width * height);
        std::vector<unsigned> row_ind;
        indInitialize(row_ind, height - kHW + 1, nHW, pHW);
        std::vector<unsigned> column_ind;
        indInitialize(column_ind, width - kHW + 1, nHW, pHW);

        // For each possible distance, precompute inter-patches distance
        for (unsigned di = 0; di <= nHW; di++) {
            for (unsigned dj = 0; dj < Ns; dj++) {
                int dk = (int) (di * width + dj) - (int) nHW;
                unsigned ddk = di * Ns + dj;

                // Process the image containing the square distance between pixels
                for (unsigned i = nHW; i < height - nHW; i++) {
                    unsigned k = i * width + nHW;
                    for (unsigned j = nHW; j < width - nHW; j++, k++)
                        diff_table[k] = (img[k + dk] - img[k]) * (img[k + dk] - img[k]);
                }

                // Compute the sum for each patches, using the method of the integral images
                unsigned dn = nHW * width + nHW;
                // 1st patch, top left corner
                float value = 0.0f;
                for (unsigned p = 0; p < kHW; p++) {
                    unsigned pq = p * width + dn;
                    for (unsigned q = 0; q < kHW; q++, pq++)
                        value += diff_table[pq];
                }
                sum_table[ddk][dn] = value;

                // 1st row, top
                for (unsigned j = nHW + 1; j < width - nHW; j++) {
                    unsigned ind = nHW * width + j - 1;
                    float sum = sum_table[ddk][ind];
                    for (unsigned p = 0; p < kHW; p++)
                        sum += diff_table[ind + p * width + kHW] - diff_table[ind + p * width];
                    sum_table[ddk][ind + 1] = sum;
                }

                // General case
                for (unsigned i = nHW + 1; i < height - nHW; i++) {
                    unsigned ind = (i - 1) * width + nHW;
                    float sum = sum_table[ddk][ind];
                    // 1st column, left
                    for (unsigned q = 0; q < kHW; q++)
                        sum += diff_table[ind + kHW * width + q] - diff_table[ind + q];
                    sum_table[ddk][ind + width] = sum;

                    // Other columns
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

        // Precompute Bloc Matching
        std::vector<std::pair<float, unsigned> > table_distance;
        table_distance.reserve(Ns * Ns);

        for(unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++) {
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
                unsigned k_r = row_ind[ind_i] * width + column_ind[ind_j];
                table_distance.clear();
                patch_table[k_r].clear();

                // Threshold distances in order to keep similar patches
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
                                        closestPower2(table_distance.size()) : NHW);

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

    void bior2dInverseProcess(std::vector<float> &group_3D_table, unsigned kHW) const {
        unsigned kHW_2 = kHW * kHW;
        unsigned N = group_3D_table.size() / kHW_2;

        // Bior process
        for (unsigned n = 0; n < N; n++)
            bior2dInverse(group_3D_table, kHW, n * kHW_2);
    }

    void bior2dProcess(std::vector<float> &bior_table_2D, std::vector<float> &img, unsigned nHW, 
        unsigned width, unsigned height, unsigned chnls, unsigned kHW, unsigned i_r, unsigned step,
        unsigned i_min, unsigned i_max) const {
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
                        bior2dForward(img, bior_table_2D, kHW, dc +
                                 (i_r + i - nHW) * width + j, width,
                                 dc_p + (i * width + j) * kHW_2);
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
                        bior2dForward(img, bior_table_2D, kHW, dc +
                                 (i + 2 * nHW + 1 - step + i_r - nHW) * width + j,
                                 width, dc_p + ((i + 2 * nHW + 1 - step)
                                 * width + j) * kHW_2);
                    }
            }
        }
    }

    void htFilteringHadamard(std::vector<float> &group_3D, std::vector<float> &tmp, unsigned nSx_r, unsigned kHard,
        unsigned chnls, std::vector<float> &sigma_table, float lambdaHard3D, std::vector<float> &weight_table) const {
        unsigned kHard_2 = kHard * kHard;
        for (unsigned c = 0; c < chnls; c++)
            weight_table[c] = 0.0f;
        float coef_norm = sqrtf((float) nSx_r);
        float coef = 1.0f / (float) nSx_r;

        // Process the Welsh-Hadamard transform on the 3rd dimension
        for (unsigned n = 0; n < kHard_2 * chnls; n++)
            hadamardTransform(group_3D, tmp, nSx_r, n * nSx_r);

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
            hadamardTransform(group_3D, tmp, nSx_r, n * nSx_r);

        for (unsigned k = 0; k < group_3D.size(); k++)
            group_3D[k] *= coef;
    }

    void wienerFilteringHadamard(std::vector<float> &group_3D_img, std::vector<float> &group_3D_est, std::vector<float> &tmp,
        unsigned nSx_r, unsigned kWien, unsigned chnls, std::vector<float> &sigma_table, std::vector<float> &weight_table) const {
        unsigned kWien_2 = kWien * kWien;
        float coef = 1.0f / (float) nSx_r;

        for (unsigned c = 0; c < chnls; c++)
            weight_table[c] = 0.0f;

        // Process the Welsh-Hadamard transform on the 3rd dimension
        for (unsigned n = 0; n < kWien_2 * chnls; n++) {
            hadamardTransform(group_3D_img, tmp, nSx_r, n * nSx_r);
            hadamardTransform(group_3D_est, tmp, nSx_r, n * nSx_r);
        }

        // Wiener Filtering
        for (unsigned c = 0; c < chnls; c++) {
            unsigned dc = c * nSx_r * kWien_2; 
            for (unsigned k = 0; k < kWien_2 * nSx_r; k++) {
                float value = group_3D_est[dc + k] * group_3D_est[dc + k] * coef;
                value /= (value + sigma_table[c] * sigma_table[c]);
                group_3D_est[k + dc] = group_3D_img[k + dc] * value * coef;
                weight_table[c] += (value*value);
            }
        }

        // Process of the Welsh-Hadamard inverse transform
        for (unsigned n = 0; n < kWien_2 * chnls; n++)
            hadamardTransform(group_3D_est, tmp, nSx_r, n * nSx_r);
    }

    void sdWeighting(std::vector<float> &group_3D, unsigned nSx_r, unsigned kHW, unsigned chnls,
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

    void bm3d1stStep(float sigma, std::vector<float> & img_noisy, std::vector<float> &img_basic, unsigned width, 
        unsigned height, unsigned chnls, unsigned nHard, unsigned kHard, unsigned NHard, unsigned pHard) const {
        float lambdaHard3D = 2.7f;                // Threshold for Hard Thresholding
        float tauMatch; // Threshold used to determinate similarity between patches
        unsigned kHard_2 = kHard * kHard;
        std::vector<float> sigma_table(chnls);    // Estimatation of sigma on each channel
        std::vector<unsigned> row_ind, column_ind;// Row, Col index step by pHard
        std::vector<float> group_3D_table(chnls * kHard_2 * NHard * column_ind.size());
        std::vector<float> wx_r_table(chnls * column_ind.size());
        std::vector<float> hadamard_tmp(NHard);
        std::vector<float> denominator(width * height * chnls, 0.0f);
        std::vector<float> numerator  (width * height * chnls, 0.0f);
        std::vector<std::vector<unsigned> > patch_table;
        std::vector<float> table_2D((2 * nHard + 1) * width * chnls * kHard_2, 0.0f); // wavelet transform
        
        mapSigma(sigma, sigma_table);
        tauMatch = sigma_table[0] < 35.0f ? 2500 : 5000;
        indInitialize(row_ind, height - kHard + 1, nHard, pHard);
        indInitialize(column_ind, width - kHard + 1, nHard, pHard);
        img_basic.resize(img_noisy.size());

        // Precompute Bloc-Matching
        precomputeBM(patch_table, img_noisy, width, height, kHard, NHard, nHard, pHard, tauMatch);

        // Loop on i_r
        for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++) {
            unsigned i_r = row_ind[ind_i];

            bior2dProcess(table_2D, img_noisy, nHard, width, height, chnls,
                kHard, i_r, pHard, row_ind[0], row_ind.back());

            wx_r_table.clear();
            group_3D_table.clear();

            // Loop on j_r
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++) {
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
                htFilteringHadamard(group_3D, hadamard_tmp, nSx_r, kHard, chnls, sigma_table,
                    lambdaHard3D, weight_table);

                sdWeighting(group_3D, nSx_r, kHard, chnls, weight_table);

                // Save the 3D group. The BIOR 2D inverse will be done after.
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++)
                        for (unsigned k = 0; k < kHard_2; k++)
                            group_3D_table.push_back(group_3D[n + k * nSx_r + c * kHard_2 * nSx_r]);

                // Save weighting
                for (unsigned c = 0; c < chnls; c++)
                    wx_r_table.push_back(weight_table[c]);

            } // End of loop on j_r
            
            bior2dInverseProcess(group_3D_table, kHard);

            // Registration of the weighted estimation
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

    void bm3d2ndStep(float sigma, std::vector<float> &img_noisy, std::vector<float> &img_basic, std::vector<float> &img_denoised,
        unsigned width, unsigned height, unsigned chnls, unsigned nWien, unsigned kWien, unsigned NWien, unsigned pWien) const {
        float tauMatch;
        unsigned kWien_2 = kWien * kWien;
        std::vector<float> sigma_table(chnls);
        std::vector<unsigned> row_ind, column_ind;
        std::vector<float> group_3D_table(chnls * kWien_2 * NWien * column_ind.size());
        std::vector<float> wx_r_table(chnls * column_ind.size());
        std::vector<float> tmp(NWien);
        std::vector<float> denominator(width * height * chnls, 0.0f);
        std::vector<float> numerator  (width * height * chnls, 0.0f);
        std::vector<std::vector<unsigned> > patch_table;
        std::vector<float> table_2D_img((2 * nWien + 1) * width * chnls * kWien_2, 0.0f);
        std::vector<float> table_2D_est((2 * nWien + 1) * width * chnls * kWien_2, 0.0f);

        mapSigma(sigma, sigma_table);
        tauMatch = (sigma_table[0] < 35.0f ? 400 : 3500);
        indInitialize(row_ind, height - kWien + 1, nWien, pWien);
        indInitialize(column_ind, width - kWien + 1, nWien, pWien);
        img_denoised.resize(img_noisy.size());

        // Precompute Bloc-Matching
        precomputeBM(patch_table, img_basic, width, height, kWien, NWien, nWien, pWien, tauMatch);

        // Loop on i_r
        for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++) {
            unsigned i_r = row_ind[ind_i];

            bior2dProcess(table_2D_img, img_noisy, nWien, width, height,
                chnls, kWien, i_r, pWien, row_ind[0], row_ind.back());
            bior2dProcess(table_2D_est, img_basic, nWien, width, height,
                chnls, kWien, i_r, pWien, row_ind[0], row_ind.back());

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
                    for (unsigned n = 0; n < nSx_r; n++)
                    {
                        const unsigned ind = patch_table[k_r][n] + (nWien - i_r) * width;
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
                wienerFilteringHadamard(group_3D_img, group_3D_est, tmp, nSx_r, kWien,
                                                chnls, sigma_table, weight_table);

                // 3D weighting using Standard Deviation
                sdWeighting(group_3D_est, nSx_r, kWien, chnls, weight_table);

                // Save the 3D group. The Bior 2D inverse will be done after.
                for (unsigned c = 0; c < chnls; c++)
                    for (unsigned n = 0; n < nSx_r; n++)
                        for (unsigned k = 0; k < kWien_2; k++)
                            group_3D_table.push_back(group_3D_est[n + k * nSx_r + c * kWien_2 * nSx_r]);

                // Save weighting
                for (unsigned c = 0; c < chnls; c++)
                    wx_r_table.push_back(weight_table[c]);

            } // End of loop on j_r

            bior2dInverseProcess(group_3D_table, kWien);

            // Registration of the weighted estimation
            unsigned dec = 0;
            for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
            {
                const unsigned j_r   = column_ind[ind_j];
                const unsigned k_r   = i_r * width + j_r;
                const unsigned nSx_r = patch_table[k_r].size();
                for (unsigned c = 0; c < chnls; c++)
                {
                    for (unsigned n = 0; n < nSx_r; n++)
                    {
                        const unsigned k = patch_table[k_r][n] + c * width * height;
                        for (unsigned p = 0; p < kWien; p++)
                            for (unsigned q = 0; q < kWien; q++)
                            {
                                const unsigned ind = k + p * width + q;
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

    // Denoise the bitmap
    Bitmap *denoise(Bitmap *noisy, Sampler *sampler) const {
        std::vector<float> img_noisy, img_basic, img_denoised;
        std::vector<float> img_sym_noisy, img_sym_basic, img_sym_denoised;
        unsigned width, height, chnls;
        unsigned h_b, w_b;   // Height and width after padding
        unsigned nHard = 16; // Half size of the search window
        unsigned nWien = 16; // Half size of the search window
        unsigned NHard = 16; // Must be a power of 2
        unsigned NWien = 32; // Must be a power of 2
        unsigned pHard = 3;
        unsigned pWien = 3;
        unsigned kHard = 16;  // Kaiser window size
        unsigned kWien = 16;  // Kaiser window size
        float coef_norm = 1.f / (sqrtf(2.f) * 128.f);
        float sigma;

        loadImage(noisy, img_noisy, &width, &height, &chnls);
        noiseEstimation(img_noisy, width, height, chnls, sampler, sigma);
        cout << "Denoised sigma: " << sigma << endl;

        img_basic.resize(img_noisy.size());
        img_denoised.resize(img_noisy.size());
        h_b   = height + 2 * nHard;
        w_b   = width  + 2 * nHard;
        for (unsigned k = 0; k < 10; k++) {
            LoD[k] *= coef_norm;
            HiR[k] *= coef_norm;
        }

        // Transformation to OPP color space
        colorSpaceTransform(img_noisy, width, height, chnls, true);

        symmetricPadding(img_noisy, img_sym_noisy, width, height, chnls, nHard);
        bm3d1stStep(sigma, img_sym_noisy, img_sym_basic, w_b, h_b, chnls, nHard, kHard, NHard, pHard);
        for (unsigned c = 0; c < chnls; c++) {
            unsigned dc_b = c * w_b * h_b + nHard * w_b + nHard;
            unsigned dc = c * width * height;
            for (unsigned i = 0; i < height; i++)
                for (unsigned j = 0; j < width; j++, dc++)
                    img_basic[dc] = img_sym_basic[dc_b + i * w_b + j];
        }

        symmetricPadding(img_basic, img_sym_basic, width, height, chnls, nHard);
        bm3d2ndStep(sigma, img_sym_noisy, img_sym_basic, img_sym_denoised, w_b, h_b, chnls, nWien, kWien, NWien, pWien); 
        for (unsigned c = 0; c < chnls; c++) {
            unsigned dc_b = c * w_b * h_b + nWien * w_b + nWien;
            unsigned dc = c * width * height;
            for (unsigned i = 0; i < height; i++)
                for (unsigned j = 0; j < width; j++, dc++)
                    img_denoised[dc] = img_sym_denoised[dc_b + i * w_b + j];
        }

        // Inverse color space transform to RGB
        colorSpaceTransform(img_denoised, width, height, chnls, false);
        colorSpaceTransform(img_noisy, width, height, chnls, false);
        colorSpaceTransform(img_basic, width, height, chnls, false);

        Bitmap *denoisedResults = storeBitmap(img_denoised, width, height);
        return denoisedResults;
    }
    std::string toString() const {
        return "BM3D denoiser\n";
    }
};

NORI_REGISTER_CLASS(BM3D_denoiser, "bm3d_denoiser");
NORI_NAMESPACE_END
