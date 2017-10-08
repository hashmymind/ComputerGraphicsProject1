#include "Application.h"
#include "qt_opengl_framework.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>

using Pixel = std::tuple<unsigned char, unsigned char, unsigned char>;

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene(void)
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	//openImage("0.png");

}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage(QString filePath)
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();
	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath)
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB(void)
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (!img_data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0; j < img_width; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0; i < 3; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k = 0; k<3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			img_data[offset_rgba + rr] = (rgb[offset_rgb + rr] >> 5) << 5;
			img_data[offset_rgba + gg] = (rgb[offset_rgb + gg] >> 5) << 5;
			img_data[offset_rgba + bb] = (rgb[offset_rgb + bb] >> 5) << 5;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned int pixelCount = img_width * img_height;
	unsigned char *rgb = this->To_RGB();
	std::vector<Pixel> pixels(pixelCount);
	//
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int pixelIdx = i*img_width + j;
			rgb[offset_rgb + rr] = (rgb[offset_rgb + rr] >> 3) << 3;
			rgb[offset_rgb + gg] = (rgb[offset_rgb + gg] >> 3) << 3;
			rgb[offset_rgb + bb] = (rgb[offset_rgb + bb] >> 3) << 3;
			pixels[pixelIdx] = std::make_tuple(rgb[offset_rgb + rr], rgb[offset_rgb + gg], rgb[offset_rgb + bb]);
		}
	}
	//
	std::map<Pixel, int> countColor;
	for (int i = 0; i < pixelCount; ++i) {
		countColor[pixels[i]]++;
	}
	std::vector<std::pair<Pixel, int> > vpCount(countColor.begin(), countColor.end());
	//std::copy(countColor.begin(), countColor.end(), back_inserter(vpCount));
	std::sort(vpCount.begin(), vpCount.end(), [](const std::pair<Pixel, int>& a, const std::pair<Pixel, int>& b) {return a.second > b.second; });
	//
	std::map<Pixel, Pixel> closest;
	unsigned int bound = 256 < vpCount.size() ? 256 : vpCount.size();
	for (auto it : countColor) {
		unsigned int minIdx = 0, minVal = 0xffffffff;
		for (int j = 0; j < bound; ++j) {
			int dist = pow(std::get<0>(it.first) - std::get<0>(vpCount[j].first), 2);
			dist += pow(std::get<1>(it.first) - std::get<1>(vpCount[j].first), 2);
			dist += pow(std::get<2>(it.first) - std::get<2>(vpCount[j].first), 2);
			if (dist < minVal) {
				minVal = dist;
				minIdx = j;
			}
			if (!dist)break;
		}
		closest[it.first] = vpCount[minIdx].first;
	}


	for (int i = 0; i < img_height; ++i)
	{
		for (int j = 0; j < img_width; ++j)
		{
			int pixelIdx = i*img_width + j;
			int offset_rgba = i*img_width * 4 + j * 4;
			Pixel closeColor = closest[pixels[pixelIdx]];
			img_data[offset_rgba + rr] = std::get<0>(closeColor);
			img_data[offset_rgba + gg] = std::get<1>(closeColor);
			img_data[offset_rgba + bb] = std::get<2>(closeColor);
			img_data[offset_rgba + aa] = WHITE;
		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}





void Application::MedianCut() {
	void;
}






//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();
	const unsigned char threshold = 127;
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			unsigned char result = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			result = result > threshold ? 255 : 0;
			img_data[offset_rgba + rr] = result;
			img_data[offset_rgba + gg] = result;
			img_data[offset_rgba + bb] = result;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();

	const unsigned char threshold = 127, randMax = 51;// 51 = 255*0.2
	int randVal;

	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			randVal = (rand() % (randMax * 2)) - randMax;
			int result = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			result += randVal;
			result = result > threshold ? 255 : 0;
			img_data[offset_rgba + rr] = result;
			img_data[offset_rgba + gg] = result;
			img_data[offset_rgba + bb] = result;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();
	int briSize = img_width*img_height;
	float threshold = 0.5;
	std::vector<float> bright(briSize);
	for (int i = 0; i < briSize; ++i) {
		bright[i] = 0.3 * rgb[i * 3 + rr] + 0.59 * rgb[i * 3 + gg] + 0.11 * rgb[i * 3 + bb];
		bright[i] /= 255.0;
	}
	for (int i = 0; i < img_height; i++)
	{
		for (int j = (i & 1) ? img_width - 1 : 0; (i & 1) ? j >= 0 : j<img_width; j += (i & 1) ? -1 : 1)// z path
																										//for(int j=0;j<img_width;++j)
		{
			int offset_rgba = i*img_width * 4 + j * 4;
			int offset_bright = i*img_width + j;
			bool result = bright[offset_bright] > threshold ? true : false;
			float error = result ? (bright[offset_bright] - 1) : bright[offset_bright];
			if (j != img_width - 1)bright[offset_bright + 1] += 0.4375*error;
			if (j != 0 && i != img_height - 1)bright[offset_bright + img_width - 1] += 0.1875*error;
			if (i != img_height - 1)bright[offset_bright + img_width] += 0.3125*error;
			if (i != img_height - 1 && j != img_width - 1)bright[offset_bright + img_width + 1] += 0.0625*error;
			img_data[offset_rgba + rr] = result ? 255 : 0;
			img_data[offset_rgba + gg] = result ? 255 : 0;
			img_data[offset_rgba + bb] = result ? 255 : 0;
			img_data[offset_rgba + aa] = WHITE;


		}
	}


	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
	float threshold;
	std::vector<float> grayList(img_width*img_height);
	unsigned int graySum = 0;
	double brightness = 0;
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_pixel = i*img_width + j;
			grayList[offset_pixel] = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			graySum += grayList[offset_pixel];
		}
	}
	brightness = graySum / (double)(img_width*img_height);
	brightness /= 255.0;
	std::vector<float> sortedGrayList(grayList);
	std::sort(sortedGrayList.begin(), sortedGrayList.end());
	threshold = sortedGrayList[(1 - brightness)*(sortedGrayList.size() - 1)] / 255.0;
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_pixel = i*img_width + j;
			int offset_rgba = i*img_width * 4 + j * 4;
			unsigned char result = (grayList[offset_pixel] / 255.0) > threshold ? 255 : 0;
			img_data[offset_rgba + rr] = result;
			img_data[offset_rgba + gg] = result;
			img_data[offset_rgba + bb] = result;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();
	float clusterthreshold[4][4] = {
		{ .7059,.3529,.5882,.2353 },
		{ .0588,.9412,.8235,.4118 },
		{ .4706,.7647,.8824,.1176 },
		{ .1765,.5294,.2941,.6471 }
	};
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width * 4 + j * 4;
			int result = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];
			result = (result / (float)255 > clusterthreshold[i % 4][j % 4]) ? 255 : 0;
			img_data[offset_rgba + rr] = result;
			img_data[offset_rgba + gg] = result;
			img_data[offset_rgba + bb] = result;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
	unsigned char colorTable[3][8] = {
		{ 0,36,73,109,146,182,219,255 },
		{ 0,36,73,109,146,182,219,255 },
		{ 0,85,170,255 }
	};
	int tmp;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = (i & 1) ? img_width - 1 : 0; (i & 1) ? j >= 0 : j<img_width; j += (i & 1) ? -1 : 1)// z path
																										//for(int j=0;j<img_width;++j)
		{
			int offset_rgb = i*img_width * 3 + j * 3, offset_othergb = 0;
			int offset_rgba = i*img_width * 4 + j * 4;
			int errors[3] = { 0x7fffffff,0x7fffffff,0x7fffffff };
			unsigned char newRGB[3];
			for (int k = 0; k < 8; ++k) {
				int errorVal = rgb[offset_rgb + rr] - colorTable[0][k];
				if (std::abs(errorVal) < std::abs(errors[0])) {
					newRGB[0] = colorTable[0][k];
					errors[0] = errorVal;
				}
				errorVal = rgb[offset_rgb + gg] - colorTable[1][k];
				if (std::abs(errorVal) < std::abs(errors[1])) {
					newRGB[1] = colorTable[1][k];
					errors[1] = errorVal;
				}
				if (k < 4) {
					errorVal = rgb[offset_rgb + bb] - colorTable[2][k];
					if (std::abs(errorVal) < std::abs(errors[2])) {
						newRGB[2] = colorTable[2][k];
						errors[2] = errorVal;
					}
				}
			}
			if (j != img_width - 1) {
				offset_othergb = i*img_width * 3 + (j + 1) * 3;
				for (int k = 0; k < 3; ++k) {
					tmp = rgb[offset_othergb + 2 - k] + 0.4375*errors[k];
					tmp = tmp > 255 ? 255 : tmp;
					tmp = tmp < 0 ? 0 : tmp;
					rgb[offset_othergb + 2 - k] = tmp;
				}
			}
			if (j != 0 && i != img_height - 1) {
				offset_othergb = (i + 1)*img_width * 3 + (j - 1) * 3;
				for (int k = 0; k < 3; ++k) {
					tmp = rgb[offset_othergb + 2 - k] + 0.1875*errors[k];
					tmp = tmp > 255 ? 255 : tmp;
					tmp = tmp < 0 ? 0 : tmp;
					rgb[offset_othergb + 2 - k] = tmp;
				}
			}
			if (i != img_height - 1) {
				offset_othergb = (i + 1)*img_width * 3 + j * 3;
				for (int k = 0; k < 3; ++k) {
					tmp = rgb[offset_othergb + 2 - k] + 0.3125*errors[k];
					tmp = tmp > 255 ? 255 : tmp;
					tmp = tmp < 0 ? 0 : tmp;
					rgb[offset_othergb + 2 - k] = tmp;
				};
			}
			if (i != img_height - 1 && j != img_width - 1) {
				offset_othergb = (i + 1)*img_width * 3 + (j + 1) * 3;
				for (int k = 0; k < 3; ++k) {
					tmp = rgb[offset_othergb + 2 - k] + 0.0625*errors[k];
					tmp = tmp > 255 ? 255 : tmp;
					tmp = tmp < 0 ? 0 : tmp;
					rgb[offset_othergb + 2 - k] = tmp;
				};
			}
			img_data[offset_rgba + rr] = newRGB[0];
			img_data[offset_rgba + gg] = newRGB[1];
			img_data[offset_rgba + bb] = newRGB[2];
			img_data[offset_rgba + aa] = WHITE;
		}
	}



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering(double filter[][5])
{
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb;
			int offset_rgba = i*img_width * 4 + j * 4;
			float newColor[3] = { 0 };
			for (int ni = i - 2; ni <= i + 2; ++ni) {
				if (ni >= img_height || ni < 0)continue;
				for (int nj = j - 2; nj <= j + 2; ++nj) {
					if (nj >= img_width || nj < 0)continue;
					offset_rgb = ni*img_width * 3 + nj * 3;
					newColor[0] += rgb[offset_rgb + rr] * filter[ni - i + 2][nj - j + 2];
					newColor[1] += rgb[offset_rgb + gg] * filter[ni - i + 2][nj - j + 2];
					newColor[2] += rgb[offset_rgb + bb] * filter[ni - i + 2][nj - j + 2];
				}
			}
			for (int k = 0; k < 3; ++k) {
				if (newColor[k] > 255)newColor[k] = 255;
				if (newColor[k] < 0) newColor[k] = 0;
			}
			img_data[offset_rgba + rr] = newColor[0];
			img_data[offset_rgba + gg] = newColor[1];
			img_data[offset_rgba + bb] = newColor[2];
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::filtering(double **filter, int n)
{
	unsigned char *rgb = this->To_RGB();
	int subt = n / 2;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb;
			int offset_rgba = i*img_width * 4 + j * 4;
			float newColor[3] = { 0 };

			for (int ni = i - subt; ni < i-subt+n; ++ni) {
				if (ni >= img_height || ni < 0)continue;
				for (int nj = j - subt; nj < j - subt + n; ++nj) {
					if (nj >= img_width || nj < 0)continue;
					offset_rgb = ni*img_width * 3 + nj * 3;
					newColor[0] += rgb[offset_rgb + rr] * filter[ni - i + subt][nj - j + subt];
					newColor[1] += rgb[offset_rgb + gg] * filter[ni - i + subt][nj - j + subt];
					newColor[2] += rgb[offset_rgb + bb] * filter[ni - i + subt][nj - j + subt];
				}
			}
			for (int k = 0; k < 3; ++k) {
				if (newColor[k] > 255)newColor[k] = 255;
				if (newColor[k] < 0) newColor[k] = 0;
			}
			if (n % 2 == 0 && i-1>0 && j-1>0) {
				for (int k = -1; k <= 0; ++k)for (int l = -1; l <= 0; ++l) {
					offset_rgba = (i+k)*img_width * 4 + (j+l) * 4;
					img_data[offset_rgba + rr] = newColor[0];
					img_data[offset_rgba + gg] = newColor[1];
					img_data[offset_rgba + bb] = newColor[2];
					img_data[offset_rgba + aa] = WHITE;
				}
			}
			else {
				img_data[offset_rgba + rr] = newColor[0];
				img_data[offset_rgba + gg] = newColor[1];
				img_data[offset_rgba + bb] = newColor[2];
				img_data[offset_rgba + aa] = WHITE;
			}
		}
	}

	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double weight[5][5] = {
		{ 1,1,1,1,1 },
		{ 1,1,1,1,1 },
		{ 1,1,1,1,1 },
		{ 1,1,1,1,1 },
		{ 1,1,1,1,1 }
	};
	for (int i = 0; i < 5; ++i)for (int j = 0; j < 5; ++j) {
		weight[i][j] /= 25.0;
	}
	filtering(weight);
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double weight[5][5] = {
		{ 1,2,3,2,1 },
		{ 2,4,6,4,2 },
		{ 3,6,9,6,3 },
		{ 2,4,6,4,2 },
		{ 1,2,3,2,1 }
	};
	for (int i = 0; i < 5; ++i)for (int j = 0; j < 5; ++j) {
		weight[i][j] /= 81.0;
	}
	filtering(weight);
	int pixelCount = img_width * img_height;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	double weight[5][5] = {
		{ 1, 4, 6, 4,1 },
		{ 4,16,24,16,4 },
		{ 6,24,36,24,6 },
		{ 4,16,24,16,4 },
		{ 1, 4, 6, 4,1 }
	};
	for (int i = 0; i < 5; ++i)for (int j = 0; j < 5; ++j) {
		weight[i][j] /= 256.0;
	}
	filtering(weight);
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N(unsigned int N)
{
	double **weight = new double*[N];

	for (int i = 0; i < N; ++i) weight[i] = new double[N]();
	double sum = 0.0;
	for (int r = N/2; r < N; ++r) {
		weight[0][r] = 1;
		for (int i = 2; i <= (N - r - 1); ++i) {
			weight[0][r] /= i;
		}
		for (int i = (r + 1); i <= N - 1; ++i) {
			weight[0][r] *= i;
		}
		weight[0][N - r - 1] = weight[0][r];
	}
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			weight[i][j] = weight[0][i] * weight[0][j];
			sum += weight[i][j];
		}
	}
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			weight[i][j] /= sum;
		}
	}
	filtering(weight, N);
	for (int i = 0; i < N; ++i) delete[] weight[i];
	delete[] weight;
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	double weight[5][5] = {
		{ -1, -4, -6, -4,-1 },
		{ -4,-16,-24,-16,-4 },
		{ -6,-24,220,-24,-6 },
		{ -4,-16,-24,-16,-4 },
		{ -1, -4, -6, -4,-1 }
	};
	for (int i = 0; i < 5; ++i)for (int j = 0; j < 5; ++j) {
		weight[i][j] /= 256.0;
	}
	filtering(weight);
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();
	double weight[5][5] = {
		{ -1, -4, -6, -4,-1 },
		{ -4,-16,-24,-16,-4 },
		{ -6,-24,476,-24,-6 },
		{ -4,-16,-24,-16,-4 },
		{ -1, -4, -6, -4,-1 }
	};
	for (int i = 0; i < 5; ++i)for (int j = 0; j < 5; ++j) {
		weight[i][j] /= 256.0;
	}
	filtering(weight);
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	unsigned char *rgb = this->To_RGB();
	img_height2 = img_height - img_height % 2;
	img_width2 = img_width - img_width % 2;
	int pixelCount = 0;
	for (int i = 0; i<img_height2; i += 2)
	{
		for (int j = 0; j<img_width2; j += 2)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba = i*img_width2 + j * 2;
			img_data[offset_rgba + rr] = rgb[offset_rgb + rr];
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg];
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb];
			++pixelCount;
		}
	}
	img_height = img_height2 / 2;
	img_width = img_width2 / 2;
	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	unsigned char *rgb = this->To_RGB();
	unsigned char *new_img_data = new unsigned char[img_height * 2 * img_width * 2 * 4];
	for (int i = 0; i<img_height; i++)
	{
		for (int j = 0; j<img_width; j++)
		{
			int offset_rgb = i*img_width * 3 + j * 3;
			int offset_rgba;
			for (int k = 0; k < 2; ++k) {
				for (int m = 0; m < 2; ++m) {
					offset_rgba = (i * 2 + k)*img_width * 8 + (j * 2 + m) * 4;
					for (int l = 0; l < 3; ++l) {
						new_img_data[offset_rgba + l] = rgb[offset_rgb + l];
					}
					new_img_data[offset_rgba + aa] = 255;
				}
			}

		}
	}
	img_width *= 2;
	img_height *= 2;
	img_data = new_img_data;
	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize(float scale)
{
	int new_width = img_width * scale, new_height = img_height * scale;
	unsigned char *new_img_data = new unsigned char[new_height * new_width * 4];
	unsigned char *rgb = this->To_RGB();
	for (int i = 0; i < new_height; i++)
	{
		for (int j = 0; j < new_width; j++)
		{
			int offset_rgb = (int)(i / scale) * 3 * img_width + (int)(j / scale) * 3;
			int offset_rgba = i*new_width * 4 + j * 4;
			new_img_data[offset_rgba + rr] = rgb[offset_rgb + rr];
			new_img_data[offset_rgba + gg] = rgb[offset_rgb + gg];
			new_img_data[offset_rgba + bb] = rgb[offset_rgb + bb];
			new_img_data[offset_rgba + aa] = WHITE;
		}
	}
	ui_instance->ui.label->setFixedHeight(new_height);
	ui_instance->ui.label->setFixedWidth(new_width);
	delete[] rgb;
	mImageDst = QImage(new_img_data, new_width, new_height, QImage::Format_ARGB32);
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate(float angleDegrees)
{
	unsigned char *rgb = this->To_RGB();
	float sint = sin(-angleDegrees * M_PI / 180.1), cost = cos(-angleDegrees * M_PI / 180.1);
	float halfH = img_height / 2.0, halfW = img_width / 2.0;
	// clean up
	memset(img_data, 0, sizeof(unsigned char) *img_width*img_height * 4);
	//
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb;
			int offset_rgba = i*img_width * 4 + j * 4;
			// convert the coordinate and rotate
			float nx = cost * (j - halfW) + sint * (halfH - i), ny = -sint * (j - halfW) + cost * (halfH - i);
			// convert back
			nx += halfW;
			ny = halfH - ny;
			// bound check
			if (nx >= img_width || nx < 0)
				continue;
			if (ny >= img_height || ny < 0)
				continue;
			//assign pixels
			offset_rgb = (int)ny *img_width * 3 + (int)nx * 3;
			img_data[offset_rgba + rr] = rgb[offset_rgb + rr];
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg];
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb];
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge(QString filePath)
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image(int tMethod)
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32);
	renew();
}

void Application::NPR_Paint_Layer(unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke(const Stroke& s)
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++)
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++)
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height))
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared)
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				}
				else if (dist_squared == radius_squared + 1)
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
	radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}


