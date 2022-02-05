#include <iostream>

#include <cmath>
#include <ctime>
#include "Image.h"
#include "basicImageManipulation.h"

#include "a9.h"

using namespace std;

void testSomeFunction() {
  // this is a good test case!
    someFunction();
};

void testBrush1() {

  Image im(1000, 1000, 3);

	vector<float> color {1.0, 0.5, 0}; // orange splats
	Image texture("Input/brush.png");

	brush(im, 500, 250, color, texture);

	im.write("Output/test_brush1.png");

}

void testBrush2() {

  Image im(1000, 1000, 3);

	vector<float> color {1.0, 1.0, 1.0}; // white
	Image texture(50, 50, 3);
  for (int n = 0; n < texture.width()*texture.height()*texture.channels(); n++) {
    texture(n) = 1;
  }

	brush(im, 500, 250, color, texture);
  brush(im, 500, 10, color, texture);
  brush(im, 500, 500, color, texture);
  brush(im, 500, 1000-25, color, texture);
  brush(im, 200, 1000-26, color, texture);

	im.write("Output/test_brush2.png");

}

void testBrush3() {

  srand(0);

  Image im(1000, 1000, 3);

	vector<float> color {1.0, 0.5, 0}; // orange splats
	Image texture("Input/brush.png");

  for (int n = 0; n < 10; n++) {
    brush(im, rand()%im.height()-50+25, rand()%im.width()-50+25, color, texture);
  }

	im.write("Output/test_brush3.png");

}

void testSingleScalePaint1() {

  srand(0);

  Image im("Input/china.png");
  // Image im("Input/DSC_8268.png");

	Image texture("Input/brush.png");

  Image out = Image(im.width(), im.height(), im.channels());

  Image importance = out + 1.0;

  singleScalePaint(im, out, importance, texture);

  out.write("Output/singleScalePaint1.png");

}

void testPainterly1() {

  srand(0);

  Image im("Input/villeperdue.png");

  Image texture("Input/brush.png");

  Image out = painterly(im, texture);

  out.write("Output/painterly1.png");

}

void testPainterly2() {

  srand(0);

  Image im("Input/castle.png");

  Image texture("Input/brush.png");

  Image out = painterly(im, texture);

  out.write("Output/painterly2.png");

}

void testTensor1() {

  // load images
  Image stata1("./Input/stata-1.png");
  Image stata2("./Input/stata-2.png");

  // compute tensors
  Image tensor1 = computeTensor(stata1);
  float maxi = tensor1.max();
  tensor1 = tensor1 / maxi;
  tensor1.write("./Output/stataTensor1.png");

  Image tensor2 = computeTensor(stata2);
  maxi = tensor2.max();
  tensor2 = tensor2 / maxi;
  tensor2.write("./Output/stataTensor2.png");

}

void testAngles1() {

  // load images
  Image stata1("./Input/angle.png");

  // compute tensors
  Image tensor1 = computeAngles(stata1);
  float maxi = tensor1.max();
  if (maxi != 0) {
    tensor1 = tensor1 / maxi;
    tensor1.write("./Output/angle.png");
  }

}

void testAngles2() {

  // load images
  Image stata1("./Input/stata-1.png");
  Image stata2("./Input/stata-2.png");

  // compute tensors
  Image tensor1 = computeAngles(stata1);
  float maxi = tensor1.max();
  if (maxi != 0) {
    tensor1 = tensor1 / maxi;
    tensor1.write("./Output/stataAngle1.png");
  }

  Image tensor2 = computeAngles(stata2);
  maxi = tensor2.max();
  if (maxi != 0) {
    tensor2 = tensor2 / maxi;
    tensor2.write("./Output/stataAngle2.png");
  }

}

void testAngles3() {

  // load images
  Image stata1("./Input/china.png");

  // compute tensors
  Image tensor1 = computeAngles(stata1);
  float maxi = tensor1.max();
  if (maxi != 0) {
    tensor1 = tensor1 / maxi;
    tensor1.write("./Output/chinaAngle.png");
  }

}

void testOrientedSingleScalePaint1() {

  srand(0);

  Image im("Input/angle.png");
  // Image im("Input/DSC_8268.png");

	Image texture("Input/longBrush2.png");

  Image out = Image(im.width(), im.height(), im.channels());

  Image importance = out + 1.0;

  singleScaleOrientedPaint(im, out, importance, texture);

  out.write("Output/orientedSingleScalePaint1.png");

}

void testOrientedSingleScalePaint2() {

  srand(0);

  Image im("Input/china.png");
  // Image im("Input/DSC_8268.png");

	Image texture("Input/longBrush2.png");

  Image out = Image(im.width(), im.height(), im.channels());

  Image importance = out + 1.0;

  singleScaleOrientedPaint(im, out, importance, texture);

  out.write("Output/orientedSingleScalePaint2.png");

}

void testOrientedPainterly1() {

  srand(0);

  Image im("Input/angle.png");

  Image texture("Input/longBrush.png");

  Image out = orientedPaint(im, texture);

  out.write("Output/orientedPainterly1.png");

}

void testOrientedPainterly2() {

  srand(0);

  Image im("Input/liz.png");

  Image texture("Input/longBrush2.png");

  Image out = orientedPaint(im, texture);

  out.write("Output/orientedPainterly2.png");

}

void testOrientedPainterly3() {

  srand(0);

  Image im("Input/china.png");

  Image texture("Input/longBrush2.png");

  Image out = orientedPaint(im, texture);

  out.write("Output/orientedPainterly3.png");

}

void testBrush() {
  testBrush1();
  testBrush2();
  testBrush3();
}

void testSingleScalePaint() {
  testSingleScalePaint1();
}

void testPainterly() {
  testPainterly1();
  testPainterly2();
}

void testTensorAngles() {
  testTensor1();
  testAngles1();
  testAngles2();
  testAngles3();
}

void testOrientedSingleScalePaint() {
  testOrientedSingleScalePaint1();
  testOrientedSingleScalePaint2();
}

void testOrientedPainterly() {
  testOrientedPainterly1();
  testOrientedPainterly2();
  testOrientedPainterly3();
}

void timingTest() {
  // test speed of the code

  vector<string>brushes = {"brush", "longBrush", "longBrush2"};

  for (string s: brushes) {

    Image texture = Image("./Input/"+s+".png");
    Image im = Image("./Input/stata.png");

    clock_t start = clock();

    Image out = orientedPaint(im, texture);

    clock_t end = clock();
    double duration = (end - start) * 1.0f / CLOCKS_PER_SEC;

    cout << "orientedPaint took: " << duration << "s" << endl;

    start = clock();

    Image out2 = painterly(im, texture);

    end = clock();

    duration = (end - start) * 1.0f / CLOCKS_PER_SEC;

    cout << "painterly took: " << duration << "s" << endl;

    out2.write("Output/stata-"+s+".png");
    out.write("Output/orientedStata-"+s+".png");

  }

}

void testBrushes() {
  // test different brushes I settled on with orientedPaint on 3 photos

  vector<string> brushes = {"1", "2", "3", "4", "brush", "longBrush", "longBrush2"};

  vector<string> targets = {"stata", "round", "archie"};

  for (string target: targets) {

    Image im = Image("./Input/"+target+".png");

    Image angle = computeAngles(im);
    angle = angle/angle.max();
    angle.write("Output/brushangle-"+target+".png");

    float sigma = 1.0;
    Image lumi = lumiChromi(im)[0];
    Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
    Image high_freq_lumi = lumi - blurred_lumi;
    Image lumi_energy = high_freq_lumi*high_freq_lumi;
    Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
    Image normalized_sharpness = sharpness/sharpness.max();

    normalized_sharpness.write("Output/brushsharp-"+target+".png");

    for (string s: brushes) {

      Image texture = Image("./Input/brushes/"+s+".png");

      Image out = orientedPaint(im, texture);
      out.write("Output/brushtest-"+target+"-"+s+".png");

    }

  }

}

void orientedAllPhotos() {
  // testing orientedPaint on multiple images and brushes
  // and outputs intermediate images as it produces results

  vector<string> brushes = {"1", "2", "3", "4", "brush", "longBrush", "longBrush2"};

  vector<string> targets = {"1", "2", "3", "4", "5", "6", "7"};

  for (string target: targets) {

    Image im = Image("./Input/photos/"+target+".png");

    Image angle = computeAngles(im);
    angle = angle/angle.max();
    angle.write("Output/brushangle-"+target+".png");

    float sigma = 1.0;
    Image lumi = lumiChromi(im)[0];
    Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
    Image high_freq_lumi = lumi - blurred_lumi;
    Image lumi_energy = high_freq_lumi*high_freq_lumi;
    Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
    Image normalized_sharpness = sharpness/sharpness.max();

    normalized_sharpness.write("Output/brushsharp-"+target+".png");

    for (string s: brushes) {

      Image texture = Image("./Input/brushes/"+s+".png");

      Image out = orientedPaint(im, texture);
      out.write("Output/custom-"+target+"-"+s+".png");

    }

  }

}

void testPaintMultiBrush() {
  // test painting multiple target images with the paintMultiBrush function
  // the target images are images I took and ones included in the pset

  vector<string> targets = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"};
  vector<string> targets_not_personal = {"liz", "archie", "maine", "china", "castle", "bp-1-2"};

  for (string target: targets) {

    Image im = Image("./Input/photos/"+target+".png");

    Image angle = computeAngles(im, 3);
    angle = angle/angle.max();
    angle.write("Output/brushangle-"+target+".png");

    Image t = computeTensor(angle);
    t.write("Output/tensor-angle-"+target+".png");

    Image out = paintMultiBrush(im);
    out.write("Output/multi-"+target+".png");

  }

  for (string target: targets_not_personal) {

    Image im = Image("./Input/"+target+".png");

    Image angle = computeAngles(im, 3);
    angle = angle/angle.max();
    angle.write("Output/brushangle-"+target+".png");

    Image t = computeTensor(angle);
    t.write("Output/tensor-angle-"+target+".png");

    Image out = paintMultiBrush(im);
    out.write("Output/multi-"+target+".png");

  }

}

void unsplashPhotos() {
  // test paintMultiBrush on the unsplash photos

  vector<string> targets = {"9"};

  for (string target: targets) {

    Image im = Image("./Input/unsplash_photos/unsplash_"+target+".png");

    Image out = paintMultiBrush(im);
    out.write("Output/unsplash-"+target+".png");

  }

}

void applePhotos() {
  // test paintMultiBrush on the apple photos

  vector<string> targets = {"1", "2", "3", "4"};

  for (string target: targets) {

    Image im = Image("./Input/apple_photos/"+target+".png");

    Image out = paintMultiBrush(im);
    out.write("Output/apple-"+target+".png");

  }

}

void testRecolor() {
  // test recoloring the a photo from random seeds, this is not very visually pleasing.

  // nice seeds: 21, 27, 28, 32, 38, 39, 41, 42 >> , 48, 47

  for (int r = 0; r < 5; r++) {

    Image im("Input/unsplash_photos/unsplash_1.png");

    vector<vector<float>> colors;

    int N = 10;

    srand(r);

    for (int n = 0; n < N; n++) {

      vector<float> color;
      color.push_back(((float) rand() / (RAND_MAX)));
      color.push_back(((float) rand() / (RAND_MAX)));
      color.push_back(((float) rand() / (RAND_MAX)));
      colors.push_back(color);

    }

    reduce_color_space(im, N, colors).write("Output/zz_reduce_cs_1_p2_seed="+to_string(r)+".png");

  }

}

void testRecolorModes() {
  // paint with recoloring using all 4 available modes:
  // RGB vs Lum and
  // with distance and without distance
  // see the definition of the k-means used in this problem in the write-up
  // we can either include the pixel separation from the mean in our distance
  // function or choose to exclude it.
  // the color difference function can also be on the chrominance or the
  // entire rgb image.


  Image im("Input/apple_photos/1.png");
  Image im2("Input/unsplash_photos/unsplash_6.png");

  int N = 4;

  for (size_t d = 0; d <= 1; d++) {
    for (size_t l = 0; l <= 1; l++) {
      vector<vector<float>> colors = get_color_samples(im, N, d, l);
      reduce_color_space(im, N, colors, d, l).write("Output/test_recolor_mode_d="+to_string(d)+"_l="+to_string(l)+"_a1.png");

      colors = get_color_samples(im2, N*3, d, l);
      reduce_color_space(im2, N*3, colors, d, l).write("Output/test_recolor_mode_d="+to_string(d)+"_l="+to_string(l)+"_u6.png");
    }
  }

}

void recolorNightDay() {
  // recolor a night image using the colors of a day image to get a
  // morning version of the night image, see figure 9

  Image im("Input/apple_photos/2.png");
  Image ref("Input/apple_photos/3.png");

  vector<vector<float>> colors;

  int N = 10;

  colors = get_color_samples(ref, N, false, false);

  vector<float> temp;

  // shifting the colors by 2 indices * you can try different indices for fun
  for (int k = 0; k < 2; k++) {

    temp = colors[0];

    for (int n = 0; n < N-1; n++) {
      colors[n] = colors[n+1];
    }

    colors[N-1] = temp;

  }

  colors[N-2] = {0.0f,0.0f,0.0f};

  // run the function
  reduce_color_space(im, N, colors, false, false).write("Output/2_recolor_w_3.png");

}

int main()
{
    // There are a couple of functions, and some extra code not put in a
    // function because it is meant for playing around while testing, i.e.
    // choosing custom colors etc.

    // Running all these functions together will take a ridiculously
    // long time since it runs a lot of different photos and brushes

    // I commeneted out the ones that produce more than 20 photos per run and aren't useful.

    Image imR = Image("./Input/china_cheat.png");
    Image imM = Image("./Input/liz_mask.png");
    Image imL = Image("./Input/liz.png");

    cout << "ee" << endl;

    // Image tensor2 = computeTensor(imR);
    Image tensor2 = addSubjectPainting(imR);
    // maxi = tensor2.max();
    // tensor2 = tensor2 / maxi;
    tensor2.write("./Output/testround.png");

    return EXIT_SUCCESS;

    // Basic testing functions
    testBrush();
    testSingleScalePaint();
    testPainterly();

    testTensorAngles();
    testOrientedSingleScalePaint();
    testOrientedPainterly();

    // timing test
    timingTest();
    // choosing brushes test
    testBrushes();

    // This tests orientedPaint on all the photos included
    // orientedAllPhotos();

    // test testPaintMultiBrush on all photos
    unsplashPhotos();
    applePhotos();
    testPaintMultiBrush();

    // tests recoloring from a random seed
    testRecolor();

    // tests recolroing using different distance function,
    // see comment inside code for more info
    // testRecolorModes();

    // testcase that produces figure 9
    recolorNightDay();

    // ### extract colors from an image ###

    int N = 5;

    Image imc("Input/photos/4.png");

    vector<vector<float>> colors;

    colors = get_color_samples(imc, N, false, false);

    for (vector<float> c : colors){

      for (float i : c) cout << std::hex << int(round(i*255))<<std::dec;
      cout << endl;

    }

    // ### some custom recoloring ###

    Image im("Input/apple_photos/2.png");

    // vector<vector<float>> colors;
    vector<string> colors_input;

    // here are some hex codes I used that I found cool
    colors_input = {"#FF521B", "#EDD382", "#F2F3AE", "#FC9E4F", "#020122"};
    colors_input = {"D4FBD5","84dcc6","a5ffd6","ffa69e","ECFFD6","ff7678","ff8284","ff8d8f"};
    reverse(colors_input.begin(), colors_input.end());

    colors_input = {"5F7544", "f2f8f2", "ffffff", "4480DA"};
    colors_input = {"5F7544", "f2f8f2", "4480DA", "ffffff"};

    // hex_colors is a function I wrote that converts hex colors into 0-1 rgb
    colors = hex_colors(colors_input);

    N = colors.size();

    // run the function
    reduce_color_space(im, N, colors, false, false).write("Output/custom_recoloring.png");

    return EXIT_SUCCESS;
}
