#ifndef A9_H_PHUDVTKB
#define A9_H_PHUDVTKB

#include <cmath>
#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include "matrix.h"
#include <iostream>
#include <random>

using namespace std;

// Write your declarations here, or extend the Makefile if you add source
// files
void someFunction();

void brush(Image &out, int y, int x, vector<float> color, Image &texture);

void singleScalePaint(Image &im, Image &out, Image &importance, Image &texture, int size = 10, int N = 1000, float noise = 0.3);

Image painterly(Image &im, Image &texture, int N=10000, int size=50, float noise=0.3);

Image computeTensor(const Image &im, float sigmaG = 1, float factorSigma = 4);

Image computeAngles(const Image &im, const int &channels=1);

vector<Image> rotateBrushes(Image &texture, int n);

void singleScaleOrientedPaint(const Image &im, Image &out, Image &importance, Image &texture, int size=50, int N=10000, float noise=0.3, int nAngles=36); // without angles passed in

void singleScaleOrientedPaintHelper(const Image &im, Image &out, Image &importance, Image &texture, Image &angle, int size=50, int N=10000, float noise=0.3, int nAngles=36);

Image orientedPaint(Image &im, Image &texture, int N = 10000, int size=50, float noise=0.3);

Image paintMultiBrush(Image &im, int N = 10000, int size=50, float noise=0.3);

Image addSubjectPainting(Image &im, int N = 10000, int size=50, float noise=0.3);

Image quantize(const Image &im, int bits);
Image multiGammaPaint(Image &im, int N = 10000, int size=50, float noise=0.3);

float compute_distance(vector<float> v1, vector<float> v2);
vector<vector<float>> hex_colors(vector<string> hex_vec);
bool color_sort( const vector<float>& c1, const vector<float>& c2 );

Eigen::ArrayXXd k_means(const Eigen::ArrayXXd &A, uint16_t k, size_t iters, bool dim5=false);
Image reduce_color_space(Image &im, int samples, vector<vector<float>> color_overwrite={}, bool use_distance=false, bool chrom=false);
vector<vector<float>> get_color_samples(Image &im, int samples, bool use_distance, bool chrom);

#endif /* end of include guard: A9_H_PHUDVTKB */
