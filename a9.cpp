#include <iostream>

#include "a9.h"

using namespace std;

// Write your implementations here, or extend the Makefile if you add source
// files
void someFunction() {
    cout << "ok, that's a function" << endl;
}

void brush(Image &out, int y, int x, vector<float> color, Image &texture) {

  if ((x > texture.width()/2)  && (x < (out.width()  - texture.width()/2))  &&
      (y > texture.height()/2) && (y < (out.height() - texture.height()/2))    ) {

    for (int z=0; z < texture.channels(); z++) {
      for (int dy=-texture.height()/2; dy < texture.height()/2; dy++) {
        for (int dx=-texture.width()/2; dx < texture.width()/2; dx++) {

          float opacity = texture(texture.width()/2+dx, texture.height()/2+dy, z);

          out(x+dx, y+dy, z) = out(x+dx, y+dy, z) * (1 - opacity) +
                               color[z] * opacity;

          // out(x+dx, y+dy, z) = 1;

    }}}

  }

}

void singleScalePaint(Image &im, Image &out, Image &importance, Image &texture, int size, int N, float noise) {

  float sf = float(size) / max(texture.height(), texture.width());
  Image stexture = scaleLin(texture, sf);

  vector<float> color(im.channels());

  int splats = 0;

  while (N > splats) {

    int x = rand() % im.width();
    int y = rand() % im.height();

    float r = float(rand()) / RAND_MAX;

    if (r < importance(x, y)) {

      for (int c = 0; c < im.channels(); c++){

        float n = float(rand()) / RAND_MAX;

        color[c] = im(x,y,c)*(1 - noise/2 + n*noise);

      }

      brush(out, y, x, color, stexture);

      splats += 1;

    }

	}

}

Image painterly(Image &im, Image &texture, int N, int size, float noise) {

  Image out(im.width(), im.height(), im.channels());

  Image ones(im.width(), im.height(), im.channels());
  ones = ones + 1;

  float sigma = 1.0;

  Image lumi = lumiChromi(im)[0];
	Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
	Image high_freq_lumi = lumi - blurred_lumi;
	Image lumi_energy = high_freq_lumi*high_freq_lumi;
	Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
	Image normalized_sharpness = sharpness/sharpness.max();

  // normalized_sharpness.write("Output/sharpened_mask.png");

  singleScalePaint(im, out, ones, texture, size, N, noise);

  if (false) {

    out.write("Output/LF.png");
    Image HF(im.width(), im.height(), im.channels());
    singleScalePaint(im, HF, normalized_sharpness, texture, size/4, N, noise);
    HF.write("Output/HF.png");

  }

	singleScalePaint(im, out, normalized_sharpness, texture, size/4, N, noise);

  return out;

}


// REUSED FROM PANORAMA
Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute xx/xy/yy Tensor of an image. (stored in that order)

  vector<Image> lc = lumiChromi(im);
  Image lum = lc[0];
  Image ch = lc[1];

  Image lum_blurred = gaussianBlur_separable(lum, sigmaG);
  Image lum_gradX = gradientX(lum_blurred);
  Image lum_gradY = gradientY(lum_blurred);

  //Channels in T are IxIx, IxIy, IyIy respectively
  Image T(im.width(), im.height(), 3);

  for (int x = 0; x < im.width(); x++){
  	for (int y = 0; y < im.height(); y++){

  		T(x,y,0) = pow(lum_gradX(x, y), 2);
  		T(x,y,1) =     lum_gradX(x, y)  *  lum_gradY(x, y);
  		T(x,y,2) = pow(lum_gradY(x, y), 2);

  	}
  }

  Image structure_tensor = gaussianBlur_separable(T, sigmaG * factorSigma);

  return structure_tensor;

}

Image computeAngles(const Image &im, const int &channels) {

  Image out(im.width(), im.height(), channels);
  Image tensor = computeTensor(im);

  for (int y = 0; y < im.height(); y++) {
    for (int x = 0; x < im.width(); x++) {

      Matrix A(2, 2);
      A << tensor(x,y,0), tensor(x,y,1),
           tensor(x,y,1), tensor(x,y,2);

      Eigen::SelfAdjointEigenSolver<Matrix> sol(A);

      int minLambda = min(abs(sol.eigenvalues()[0]), abs(sol.eigenvalues()[1])) == abs(sol.eigenvalues()[1]);

      Vec2f v = sol.eigenvectors().col(minLambda);

      float angle = atan2(v[1], v[0]);

      while (angle < 0)      angle += 2*M_PI;
      while (angle > 2*M_PI) angle -= 2*M_PI;


      for (int z = 0; z < channels; z++) {
        out(x, y, z) = angle;
      }

  }}

  return out;

}

vector<Image> rotateBrushes(Image &texture, int n) {
  //  takes as input a single texture and returns a list of n textures rotated by i2π/n

	vector<Image> out;

	for (int i = 0; i < n; i++){

		float theta = -M_PI*i/n;
		out.push_back(rotate(texture, theta));

	}

	return out;

}

void singleScaleOrientedPaintHelper(const Image &im, Image &out, Image &importance, Image &texture, Image &angle, int size, int N, float noise, int nAngles) {

  float sf = float(size) / max(texture.height(), texture.width());
  Image stexture = scaleLin(texture, sf);

  vector<float> color(im.channels());

  int splats = 0;

  vector<Image> rotated_brushes = rotateBrushes(stexture, nAngles);

  while (N > splats) {

    int x = rand() % im.width();
    int y = rand() % im.height();

    float r = float(rand()) / RAND_MAX;

    if (r < importance(x, y)) {

      int theta = int(round( angle(x, y)*nAngles / (M_PI) )) % nAngles;

      for (int c = 0; c < im.channels(); c++){

        float n = float(rand()) / RAND_MAX;

        color[c] = im(x,y,c)*(1 - noise/2 + n*noise);

      }

      brush(out, y, x, color, rotated_brushes[theta]);

      splats += 1;

    }

	}

}

void singleScaleOrientedPaint(const Image &im, Image &out, Image &importance, Image &texture, int size, int N, float noise, int nAngles) {

  Image angle = computeAngles(im);

  singleScaleOrientedPaintHelper(im, out, importance, texture, angle, size, N, noise, nAngles);

}

Image orientedPaint(Image &im, Image &texture, int N, int size, float noise) {

  Image out(im.width(), im.height(), im.channels());

  Image ones(im.width(), im.height(), im.channels());
  ones = ones + 1;

  float sigma = 1.0;

  Image lumi = lumiChromi(im)[0];
	Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
	Image high_freq_lumi = lumi - blurred_lumi;
	Image lumi_energy = high_freq_lumi*high_freq_lumi;
	Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
	Image normalized_sharpness = sharpness/sharpness.max();

  // normalized_sharpness.write("Output/sharpened_mask.png");

  singleScaleOrientedPaint(im, out, ones, texture, size, N, noise);

  if (false) {
    out.write("Output/LF.png");
    Image HF(im.width(), im.height(), im.channels());
    singleScaleOrientedPaint(im, HF, normalized_sharpness, texture, size/4, N, noise);
    HF.write("Output/HF.png");
  }

	singleScaleOrientedPaint(im, out, normalized_sharpness, texture, size/4, N, noise);

  return out;

}

Image paintMultiBrush(Image &im, int N, int size, float noise) {

  Image texture("./Input/brushes/1.png");
  Image texture2("./Input/brushes/longBrush.png");
  Image texture_alt("./Input/brushes/3.png");

  Image out(im.width(), im.height(), im.channels());

  Image ones(im.width(), im.height(), im.channels());
  ones = ones + 1;

  float sigma = 1.0;

  Image lumi = lumiChromi(im)[0];
	Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
	Image high_freq_lumi = lumi - blurred_lumi;
	Image lumi_energy = high_freq_lumi*high_freq_lumi;
	Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
	Image normalized_sharpness = sharpness/sharpness.max();

  Image angles = computeAngles(im, 3);

  Image tensor = computeTensor(im);

  Image delta_angle = computeTensor(angles);

  singleScaleOrientedPaintHelper(im, out, ones, texture, angles, size, N, noise);

    // REMOVE THIS
    out.write("Output/bg.png");

    Image out2(im.width(), im.height(), im.channels());
    Image out3(im.width(), im.height(), im.channels());

    singleScaleOrientedPaintHelper(im, out2, delta_angle, texture_alt, angles, size/2, N, noise);

    out2.write("Output/delta_angle.png");

  	singleScaleOrientedPaintHelper(im, out3, tensor, texture2, angles, size/2, N, noise);

    out3.write("Output/tensor.png");
    // REMOVE THIS

  singleScaleOrientedPaintHelper(im, out, delta_angle, texture_alt, angles, size/2, N, noise);

	singleScaleOrientedPaintHelper(im, out, tensor, texture2, angles, size/2, N, noise);

  return out;

}

Image addSubjectPainting(Image &im, int N, int size, float noise) {

  Image texture("./Input/brushes/1.png");
  Image texture2("./Input/brushes/longBrush.png");
  Image texture_alt("./Input/brushes/3.png");

  Image out(im.width(), im.height(), im.channels());

  Image ones(im.width(), im.height(), im.channels());
  ones = ones + 1;

  Image angles = computeAngles(im, 3);

  Image tensor = computeTensor(im);

  Image delta_angle = computeTensor(angles);

  singleScaleOrientedPaintHelper(im, out, ones, texture, angles, size, N, noise);

  singleScaleOrientedPaintHelper(im, out, delta_angle, texture_alt, angles, size/2, N, noise);

	singleScaleOrientedPaintHelper(im, out, tensor, texture2, angles, size/2, N, noise);

  Image subject("./Input/liz.png");
  Image subject_mask("./Input/liz_mask.png");

  Image subject_angles = computeAngles(subject, 3);

  Image subject_tensor = computeTensor(subject);

  Image subject_delta_angle = computeTensor(subject_angles);

  for (int n = 0; n < subject_mask.width()*subject_mask.height()*subject_mask.channels(); n++) {
    subject_tensor(n) = subject_tensor(n)*subject_mask(n)/2;
    subject_delta_angle(n) = subject_delta_angle(n)*subject_mask(n)/2;
  }

  int dN = float(N)/(im.width()*im.height())*subject.width()*subject.height();

  singleScaleOrientedPaintHelper(subject, out, subject_mask, texture, angles, size, dN, noise);

  singleScaleOrientedPaintHelper(subject, out, subject_delta_angle, texture_alt, angles, size/2, dN, noise);

	singleScaleOrientedPaintHelper(subject, out, subject_tensor, texture2, angles, size/2, dN, noise);

  return out;

}

// Quantizes the image to 2^bits levels and scales back to 0~1
Image quantize(const Image &im, int bits) {
  // // --------- HANDOUT  PS01 ------------------------------
  // Image output(im.width(), im.height(), im.channels());
  // Quantizes the image to 2^bits levels
  // return output;

  Image res(im.width(), im.height(), im.channels());

  for (int i = 0; i < im.width() * im.height() * im.channels(); i++) {

	res(i) = round(pow(2, bits) * im(i)) / pow(2, bits);

  }

  return res;
}

Image multiGammaPaint(Image &im, int N, int size, float noise) {

  Image texture("./Input/brushes/1.png");
  Image texture2("./Input/brushes/longBrush.png");
  Image texture_alt("./Input/brushes/3.png");

  Image out(im.width(), im.height(), im.channels());

  Image ones(im.width(), im.height(), im.channels());
  ones = ones + 1;

  Image quantized = quantize(gamma_code(im, 0.25), 4);

  Image angles = computeAngles(im, 3);

  Image tensor = computeTensor(im);

  for (int y = 0; y < tensor.height(); y++) {
    for (int x = 0; x < tensor.width(); x++) {
      // normalize the x^2 and y^2 tensors
      tensor(x, y, 0) = tensor(x, y, 0)/2 + tensor(x, y, 2)/2;
      // tensor(x, y, 1) = tensor(x, y, 0)/2 + tensor(x, y, 2)/2;
      // tensor(x, y, 2) = tensor(x, y, 0)/2 + tensor(x, y, 2)/2;
    }
  }

  for (int y = 3; y < tensor.height()-3; y++) {
    for (int x = 3; x < tensor.width()-3; x++) {
      // filter out the noise in the tensor
      float sum = 0;
      for (int dy = -2; dy <= 2; dy++) {
        for (int dx = -2; dx <= 2; dx++) {
          sum += tensor(x-dx, y-dy, 0);
      }}

      // cout << sum << ",";

      if (sum > 0.15) { // 0.006 * 25 pixels
        tensor(x, y, 0) = tensor(x, y, 0);
        tensor(x, y, 1) = tensor(x, y, 0);
        tensor(x, y, 2) = tensor(x, y, 0);
      } else {
        tensor(x, y, 0) = 0;
        tensor(x, y, 1) = 0;
        tensor(x, y, 2) = 0;
      }

    }
  }

  Image delta_angle = computeTensor(angles);

  Image mean = im/2 + quantized/2;

  singleScaleOrientedPaintHelper(mean, out, ones, texture, angles, size, N, noise/3);

  singleScaleOrientedPaintHelper(mean, out, delta_angle, texture_alt, angles, size/2, N, noise/3);

	singleScaleOrientedPaintHelper(mean, out, tensor, texture2, angles, size/2, N, noise);

  return out;

}

float compute_distance(vector<float> v1, vector<float> v2) {
  // return pow( pow(v1[0]-v2[0], 2) + pow(v1[1]-v2[1], 2) + pow(v1[2]-v2[2], 2), 0.5);

  float d = 0;

  for (int n = 0; n < v1.size(); n++) {

    d += pow(v1[n]-v2[n], 2);

  }

  return pow(d, 0.5);

}

vector<vector<float>> hex_colors(vector<string> hex_vec) {

  vector<vector<float>> colors;

  for (int n = 0; n < hex_vec.size(); n++) {

      int r, g, b;

      if (hex_vec[n].size() == 7) {
        sscanf(hex_vec[n].c_str(), "#%02x%02x%02x", &r, &g, &b);
      } else {
        sscanf(hex_vec[n].c_str(), "%02x%02x%02x", &r, &g, &b);
      }

      colors.push_back( { ((float)(r)/255.0f),
                          ((float)(g)/255.0f),
                          ((float)(b)/255.0f)} );

  }

  return colors;

}

bool color_sort( const vector<float>& c1, const vector<float>& c2 ) {
   return compute_distance(c1, {0,0,0}) < compute_distance(c2, {0,0,0});
}

Eigen::ArrayXXf k_means(const Eigen::ArrayXXf &A, uint16_t k, size_t iters, bool dim5) {

  static std::random_device seed;
  static std::mt19937 rand_gen(seed());
  std::uniform_int_distribution<size_t> indxs(0, A.rows() - 1);

  const Eigen::ArrayXXf A_r = A.col(0).rowwise().replicate(k);
  const Eigen::ArrayXXf A_g = A.col(1).rowwise().replicate(k);
  const Eigen::ArrayXXf A_b = A.col(2).rowwise().replicate(k);
  const Eigen::ArrayXXf A_x = A.col(3).rowwise().replicate(k);
  const Eigen::ArrayXXf A_y = A.col(4).rowwise().replicate(k);

  Eigen::ArrayXXf means(k, A.cols());

  for (size_t cluster = 0; cluster < k; cluster++) {
    for (int n = 0; n < A.cols(); n++) {
      means(cluster, n) = A(indxs(rand_gen), n);
    }
  }

  for (size_t n = 0; n < iters; n++) {

    Eigen::ArrayXXf sum = Eigen::ArrayXXf::Zero(k, A.cols());
    Eigen::ArrayXf counts = Eigen::ArrayXf::Ones(k);

    Eigen::ArrayXXf distances;

    if (dim5==true) {
      distances = (A_r.rowwise() - means.col(0).transpose()).square() +
                  (A_g.rowwise() - means.col(1).transpose()).square() +
                  (A_b.rowwise() - means.col(2).transpose()).square() +
                  (A_x.rowwise() - means.col(3).transpose()).square() +
                  (A_y.rowwise() - means.col(4).transpose()).square();
    } else {
      distances = (A_r.rowwise() - means.col(0).transpose()).square() +
                  (A_g.rowwise() - means.col(1).transpose()).square() +
                  (A_b.rowwise() - means.col(2).transpose()).square();
    }

    for (size_t idx = 0; idx < A.rows(); idx++) {

      Eigen::ArrayXf::Index argmin;
      distances.row(idx).minCoeff(&argmin);
      sum.row(argmin) += A.row(idx).array();
      counts(argmin) += 1;

    }

    means = sum.colwise() / counts;

  }

  return means;
  
}

Image reduce_color_space(Image &im, int samples, vector<vector<float>> color_overwrite, bool use_distance, bool chrom) {

  Image out(im.width(), im.height(), im.channels());

  Image target(im.width(), im.height(), im.channels());

  vector<Image> lc = lumiChromi(im);

  if (chrom) {
    target = lc[1];
  } else {
    target = im;
  }

  if (color_overwrite.size() != 0) {
    if (color_overwrite.size() != samples) {
      cout << "NUMBER OF COLORS OVERWRITTEN INCORRECT" << endl;
      return out;
    }
  }

  Eigen::ArrayXXf A = Eigen::ArrayXXf::Zero(im.width() * im.height(), 5);

  for (int y = 0; y < im.height(); y++) {
    for (int x = 0; x < im.width(); x++) {

      int d = x + y * im.width();

      A.row(d) << im(x, y, 0), im(x, y, 1), im(x, y, 2), x/min(im.width(), im.height()), y/min(im.width(), im.height());

    }
  }

  // Perform the k_means algorithm, the bulk of this function
  Eigen::ArrayXXf means = k_means(A, samples, samples*5, use_distance);
  //

  vector<vector<float>> colors;

  for (int n = 0; n < means.rows(); n++) {

    vector<float> color = {means(n, 0), means(n, 1), means(n, 2)};

    colors.push_back(color);

  }

  // Sorting for a more determinstic color replacement
  sort(colors.begin(), colors.end(), color_sort);

  for (int y = 0; y < im.height(); y++) {
  for (int x = 0; x < im.width(); x++) {

      float d_min = 5.0f;
      int idx_d = 0;
      vector<float> vec{im(x, y, 0), im(x, y, 1), im(x, y, 2)};

      for (int i = 0; i < colors.size(); i++) {

        float dp = compute_distance(vec, colors[i]);

        if (dp <= d_min) {
          d_min = dp;
          idx_d = i;
        }

      }

      if (color_overwrite.size() == 0) {

        out(x, y, 0) = colors[idx_d][0];
        out(x, y, 1) = colors[idx_d][1];
        out(x, y, 2) = colors[idx_d][2];

      } else {

        out(x, y, 0) = color_overwrite[idx_d][0];
        out(x, y, 1) = color_overwrite[idx_d][1];
        out(x, y, 2) = color_overwrite[idx_d][2];

      }

  }}

  out = paintMultiBrush(out);

  return out;

}

vector<vector<float>> get_color_samples(Image &im, int samples, bool use_distance, bool chrom) {

    Image target(im.width(), im.height(), im.channels());

    vector<Image> lc = lumiChromi(im);

    if (chrom) {
      target = lc[1];
    } else {
      target = im;
    }

    Eigen::ArrayXXf A = Eigen::ArrayXXf::Zero(im.width() * im.height(), 5);

    for (int y = 0; y < im.height(); y++) {
      for (int x = 0; x < im.width(); x++) {

        int d = x + y * im.width();

        A.row(d) << im(x, y, 0), im(x, y, 1), im(x, y, 2), x/min(im.width(), im.height()), y/min(im.width(), im.height());

      }
    }

    // Perform the k_means algorithm, the bulk of this function
    Eigen::ArrayXXf means = k_means(A, samples, samples*5, use_distance);
    //

    vector<vector<float>> colors;

    for (int n = 0; n < means.rows(); n++) {

      vector<float> color = {means(n, 0), means(n, 1), means(n, 2)};

      colors.push_back(color);

    }

    // Sorting for a more determinstic color replacement
    sort(colors.begin(), colors.end(), color_sort);

    return colors;

}
