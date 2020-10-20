#include "Canny.h"
#include<minmax.h>
#include<iostream>
using namespace std;

#define pi 3.1415926

float filter[7][7] = {  {0.0001,0.0010,0.0021,0,-0.0021,-0.0010,-0.0001},
						{0.0014,0.0117,0.0261,0,-0.0261,-0.0117,-0.0014},
						{0.0064,0.0523,0.1171,0,-0.1171,-0.0523,-0.0064},
						{0.0106,0.0862,0.1931,0,-0.1931,-0.0862,-0.0106},
						{0.0064,0.0523,0.1171,0,-0.1171,-0.0523,-0.0064},
						{0.0014,0.0117,0.0261,0,-0.0261,-0.0117,-0.0014},
						{0.0001,0.0010,0.0021,0,-0.0021,-0.0010,-0.0001}};


CannyBase::CannyBase(unsigned short* srcImage, int width, int height, int size)
	:imageWidth(width),
	imageHeight(height),
	guassKernelSize(size){
	setSrcImage(srcImage);
};


void CannyBase::setSrcImage(unsigned short* srcImage) {
	if (imageWidth == 0 || imageHeight == 0) {
		cout << "���ȳ�ʼ��ͼ��ĳ���" << endl;
		return;
	}
	if (srcImage == NULL) {
		cout << "ͼ��Ϊ�գ�" << endl;
		return;
	}
	m_srcImage.resize(imageHeight, vector<float>(imageWidth, 0));
	for (int i = 0; i < imageHeight; ++i) {
		for (int j = 0; j < imageWidth; ++j) {
			m_srcImage[i][j] = srcImage[i * imageWidth + j];
		}
	}
}

void CannyBase::Gaussianconv() {
	if (m_srcImage.empty()) {
		cout << "m_srcImage empty !" << endl;
		return;
	}

	int centerX = guassKernelSize / 2; // 0, 1, 2, 3, 4
	int centerY = guassKernelSize / 2;

	m_imageGauss.assign(m_srcImage.begin(), m_srcImage.end());

	vector<vector<float>> GaussianKernel(guassKernelSize, vector<float>(guassKernelSize, 0));
	GetGaussianKernel(GaussianKernel);

	for (int i = centerY; i < imageHeight - centerY; i++)
	{
		//�д���(ȥ����Ե����)
		for (int j = centerX; j < imageWidth - centerX; j++)
		{
			float value = 0.0;
			//����
			for (int k = -centerY; k <= centerY; k++)
			{
				for (int l = -centerX; l <= centerX; l++)
				{
					//�����Ȩƽ��,��������ֵ
					value += m_srcImage[i][j] * GaussianKernel[k + centerY][l + centerX];
				}
			}
			if (value > 65535)
			{
				m_imageGauss[i][j] = 65535;
			}
			else {
				m_imageGauss[i][j] = value;
			}
		}
	}
}

void CannyBase::GetGaussianKernel(vector<vector<float>>& Kernel, const float sigmma) {
	int centerX = int(guassKernelSize / 2); // 0, 1, 2, 3, 4
	int centerY = int(guassKernelSize / 2);

	// ��size x size����ֵ
	double sum_value = 0.0;

	for (int i = 0; i < guassKernelSize; i++)
	{
		for (int j = 0; j < guassKernelSize; j++)
		{
			Kernel[i][j] = 1 / (2 * pi * pow(sigmma, 2)) *
				exp(-(pow((double)(i - centerX), 2) + pow((double)(j - centerY), 2)) / (2 * pow(sigmma, 2)));
			sum_value += Kernel[i][j];
		}
	}

	// gaussain���ڵ�ֵ��Ϊ1
	for (int i = 0; i < guassKernelSize; i++)
	{
		for (int j = 0; j < guassKernelSize; j++)
		{
			Kernel[i][j] = Kernel[i][j] / sum_value;
		}
	}
}

void CannyBase::SobelGradDirction() {
	if (m_imageGauss.empty()) {
		cout << "GuassImage empty ��" << endl;
		return;
	}

	m_imageSobelX.resize(imageHeight, vector<float>(imageWidth, 0));
	m_imageSobelY.resize(imageHeight, vector<float>(imageWidth, 0));
	m_SobelAmpXY.resize(imageHeight, vector<float>(imageWidth, 0));
	m_pointDrection.resize(imageHeight, vector<float>(imageWidth, 0));

	for (int i = 1; i < imageHeight - 1; ++i)
	{
		for (int j = 1; j < imageWidth - 1; ++j)
		{
			double gradX = m_imageGauss[i - 1][j + 1] + m_imageGauss[i][j + 1] * 2 + m_imageGauss[i + 1][j + 1]
				- m_imageGauss[i - 1][j - 1] - m_imageGauss[i][j - 1] * 2 - m_imageGauss[i + 1][j - 1];
			double gradY = m_imageGauss[i + 1][j - 1] + m_imageGauss[i + 1][j] * 2 + m_imageGauss[i + 1][j + 1]
				- m_imageGauss[i - 1][j - 1] - m_imageGauss[i - 1][j] * 2 - m_imageGauss[i - 1][j + 1];
			m_imageSobelY[i][j] = abs(gradY) > 65535 ? 65535 : abs(gradY);
			m_imageSobelX[i][j] = abs(gradX) > 65535 ? 65535 : abs(gradX);
			if (gradX == 0)
			{
				gradX = 0.00000000000000001;
			}
			m_pointDrection[i][j] = (atan(gradY / gradX) * 180) / pi;
			if (m_pointDrection[i][j] < 0)
				m_pointDrection[i][j] += 180;
		}
	}
	SobelAmplitude();
}

void CannyBase::SobelAmplitude()
{
	for (int i = 0; i < imageHeight; i++)
	{
		for (int j = 0; j < imageWidth; j++)
		{
			double ampXY = sqrt(m_imageSobelX[i][j] * m_imageSobelX[i][j] + m_imageSobelY[i][j] * m_imageSobelY[i][j]);
			m_SobelAmpXY[i][j] = ampXY > 65535 ? 65535 : ampXY;
		}
	}
}


//******************�Ǽ���ֵ����*************************
//��һ������imageInput�����Sobel�ݶ�ͼ��
//�ڶ�������imageOutPut������ľֲ�����ֵ����ͼ��
//���������pointDrection��ͼ����ÿ������ݶȷ�������ָ��
//*************************************************************
void CannyBase::LocalMaxValue() {
	if (m_SobelAmpXY.empty()) {
		cout << "m_SobelAmpXY empty ��" << endl;
		return;
	}

	float g1 = 0, g2 = 0, g3 = 0, g4 = 0;                            //���ڽ��в�ֵ���õ������ص�����ֵ   
	double dTmp1 = 0.0, dTmp2 = 0.0;                           //�������������ص��ֵ�õ��ĻҶ����� 
	double dWeight = 0.0;                                    //��ֵ��Ȩ��  
	m_nonMaximum.assign(m_SobelAmpXY.begin(), m_SobelAmpXY.end());

	//�趨�ĽǶ���ֵΪ45��
	for (int i = 0; i < imageHeight; i++)
	{
		for (int j = 0; j < imageWidth; j++)
		{
			if (m_SobelAmpXY[i][j] == 0) {
				continue;
			}
			else {
				//145~180
				if (m_pointDrection[i][j] >= 135 && m_pointDrection[i][j]<180)
				{
					g1 = m_SobelAmpXY[i - 1][j - 1];
					g2 = m_SobelAmpXY[i][j - 1];
					g3 = m_SobelAmpXY[i][j + 1];
					g4 = m_SobelAmpXY[i + 1][j + 1];
					dWeight = m_imageSobelY[i][j] / m_imageSobelX[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
				//0~45
				else if (m_pointDrection[i][j] >= 0 && m_pointDrection[i][j]<45)
				{
					g1 = m_SobelAmpXY[i - 1][j + 1];
					g2 = m_SobelAmpXY[i][j + 1];
					g3 = m_SobelAmpXY[i][j - 1];
					g4 = m_SobelAmpXY[i + 1][j - 1];
					dWeight = m_imageSobelY[i][j] / m_imageSobelX[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
				//45~90
				else if (m_pointDrection[i][j] >= 45 && m_pointDrection[i][j] < 90)
				{
					g1 = m_SobelAmpXY[i - 1][j + 1];
					g2 = m_SobelAmpXY[i - 1][j];
					g3 = m_SobelAmpXY[i + 1][j];
					g4 = m_SobelAmpXY[i + 1][j - 1];
					dWeight = m_imageSobelX[i][j] / m_imageSobelY[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
				else if (m_pointDrection[i][j] >= 90 && m_pointDrection[i][j] < 135)
				{
					g1 = m_SobelAmpXY[i - 1][j - 1];
					g2 = m_SobelAmpXY[i - 1][j];
					g3 = m_SobelAmpXY[i + 1][j];
					g4 = m_SobelAmpXY[i + 1][j + 1];
					dWeight = m_imageSobelX[i][j] / m_imageSobelY[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
			}
			if (m_nonMaximum[i][j] <= dTmp1 || m_nonMaximum[i][j] <= dTmp2)
				m_nonMaximum[i][j] = 0;
		}
	}
}

//˫��ֵ��ѡȡ
void CannyBase::SetThreshold() {
	int nHist[65536]{ 0 };//ֱ��ͼ
	int nEdgeNum = 0;//���б�Ե�����Ŀ
	int nMaxMag = 0;//����ݶȵķ�ֵ

	for (int i = 0; i < imageHeight; ++i) {
		for (int j = 0; j < imageWidth; ++j) {
			int nindex = m_nonMaximum[i][j];
			if (nindex > 0) {
				nHist[nindex]++;
				nEdgeNum++;
				nMaxMag = max(nindex, nMaxMag);
			}
		}
	}

	//����������ֵ ע�����ݶȵ���ֵ
	double dRateHigh = 0.5;
	double dRateLow = 0.5;
	int nHightcount = (int)(dRateHigh*nEdgeNum + 0.5);
	int count = 0;
	nEdgeNum = nHist[1];
	while ((nEdgeNum <= nHightcount) && (count < nMaxMag - 1)) {
		count++;
		nEdgeNum += nHist[count];
	}
	m_highThreshold = double(count);
	m_lowThreshold = m_highThreshold*dRateLow;
}


void CannyBase::DoubleThreshold()
{
	if (m_nonMaximum.empty()) {
		cout << "m_nonMaximum empty ��" << endl;
		return;
	}
	m_imageDouble.assign(m_nonMaximum.begin(), m_nonMaximum.end());

	//��ø���ֵ�͵���ֵ  ����ֵΪ70%������ֵΪ50%����ֵ
	SetThreshold();

	for (int i = 0; i < imageHeight; i++){
		for (int j = 0; j < imageWidth; j++){
			if (m_nonMaximum[i][j] > m_highThreshold){
				m_imageDouble[i][j] = 65535;

			}
			else if (m_nonMaximum[i][j] < m_lowThreshold){
				m_imageDouble[i][j] = 0;
			}
		}
	}
}

void CannyBase::DoubleThresholdLink()
{
	if (m_imageDouble.empty()) {
		cout << "m_imageDouble empty !" << endl;
		return;
	}

	m_imageLink.assign(m_imageDouble.begin(), m_imageDouble.end());
	vector<pair<int, int>> flag;//������¼����˫��ֵ֮��ĵ�
	int change = 0;//��¼�仯�ĵ�ĸ���
	for (int i = 0; i <imageHeight; i++){
		for (int j = 0; j < imageWidth; j++){
			if (m_imageLink[i][j] > m_lowThreshold&&m_imageLink[i][j] <= m_highThreshold){
				if (m_imageLink[i - 1][j - 1] == 65535 || m_imageLink[i - 1][j] == 65535 || m_imageLink[i - 1][j + 1] == 65535 ||
					m_imageLink[i][j - 1] == 65535 || m_imageLink[i][j + 1] == 65535 || m_imageLink[i + 1][j - 1] == 65535 ||
					m_imageLink[i + 1][j] == 65535 || m_imageLink[i + 1][j + 1] == 65535){
					m_imageLink[i][j] = 65535;
					change++;
				}
				else{
					flag.push_back(make_pair(i, j));
				}
			}
		}
	}

	while (change > 0 && !flag.empty()) {
		change = 0;
		vector<pair<int, int>> tempflag;
		for (auto it : flag)
		{
			int i = it.first;
			int j = it.second;
			if (m_imageLink[i - 1][j - 1] == 65535 || m_imageLink[i - 1][j] == 65535 || m_imageLink[i - 1][j + 1] == 65535 ||
				m_imageLink[i][j - 1] == 65535 || m_imageLink[i][j + 1] == 65535 || m_imageLink[i + 1][j - 1] == 65535 ||
				m_imageLink[i + 1][j] == 65535 || m_imageLink[i + 1][j + 1] == 65535)
			{
				m_imageLink[i][j] = 65535;
				change++;
			}
			else {
				tempflag.push_back(it);
			}
		}
		flag = tempflag;
	}
	for (auto it : flag) {
		int i = it.first;
		int j = it.second;
		m_imageLink[i][j] = 0;
	}
}

void CannyBase::CannyEdgeDetect() {
	if (m_srcImage.empty()) {
		cout << "m_srcImage empty !" << endl;
		return;
	}
	Gaussianconv();
	SobelGradDirction();
	LocalMaxValue();
	DoubleThreshold();
	DoubleThresholdLink();

}













MyCanny::MyCanny(unsigned short* srcImage, int width, int height, int size)
	:CannyBase(srcImage,width,height,size) {

};

//void MyCanny::GetGaussianKernel(vector<vector<float>>& Kernel, const int KernelSize , const float sigmma ) {
//	int centerX = int(KernelSize / 2); // 0, 1, 2, 3, 4
//	int centerY = int(KernelSize / 2);
//
//	// ��size x size����ֵ
//	double sum_value = 0.0;
//
//	for (int i = 0; i < KernelSize; i++)
//	{
//		for (int j = 0; j < KernelSize; j++)
//		{
//			Kernel [i][j] = 1 / (2 * pi * pow(sigmma, 2)) *
//				exp(-(pow((double)(i - centerX), 2) + pow((double)(j - centerY), 2)) / (2 * pow(sigmma, 2)));
//			sum_value += Kernel[i][j];
//		}
//	}
//
//	// gaussain���ڵ�ֵ��Ϊ1
//	for (int i = 0; i < KernelSize; i++)
//	{
//		for (int j = 0; j < KernelSize; j++)
//		{
//			Kernel[i][j] = Kernel[i][j] / sum_value;
//		}
//	}
//}
//
//void MyCanny::Gaussianconv(unsigned short* pImageIn, int width, int height, const int size) {
//	int centerX = size / 2; // 0, 1, 2, 3, 4
//	int centerY = size / 2;
//	//����ԭʼͼ��
//	m_srcImage.resize(height, vector<float>(width, 0));
//	for (int i = 0; i < height; ++i) {
//		for (int j = 0; j < width; ++j) {
//			m_srcImage[i][j] = pImageIn[i * width + j];
//		}
//	}
//
//	m_imageGauss.assign(m_srcImage.begin(),m_srcImage.end());
//
//	vector<vector<float>> GaussianKernel(size, vector<float>(size, 0));
//	GetGaussianKernel(GaussianKernel, size);
//
//	for (int i = centerY; i < height - centerY; i++)
//	{
//		//�д���(ȥ����Ե����)
//		for (int j = centerX; j < width - centerX; j++)
//		{
//			float value = 0.0;
//			//����
//			for (int k = -centerY; k <= centerY; k++)
//			{
//				for (int l = -centerX; l <= centerX; l++)
//				{
//					//�����Ȩƽ��,��������ֵ
//					value += pImageIn[width * i + j + l] * GaussianKernel[k +centerY][l+centerX];
//				}
//			}
//			if (value > 65535)
//			{
//				m_imageGauss[i][j] = 65535;
//			}
//			else {
//				m_imageGauss[i][j] = value;
//			}
//		}
//	}
//}


//�洢�ݶ�Ĥ�����ݶȽ�
void MyCanny::SobelGradDirction()
{
	if (m_imageGauss.empty()) {
		cout << "MyCanny GuassImage empty ��" << endl;
		return;
	}

	m_imageSobelX.resize(imageHeight, vector<float>(imageWidth, 0));
	m_imageSobelY.resize(imageHeight, vector<float>(imageWidth, 0));
	m_SobelAmpXY.resize(imageHeight, vector<float>(imageWidth, 0));
	m_pointDrection.resize(imageHeight, vector<float>(imageWidth, 0));

	for (int i = m_filterCenter; i < imageHeight - m_filterCenter; ++i)
	{
		for (int j = m_filterCenter; j < imageWidth - m_filterCenter; ++j)
		{
			double gradX = 0;
			double gradY = 0;
			for (int m = -m_filterCenter; m <= m_filterCenter; ++m) {
				for (int n = -m_filterCenter; n <= m_filterCenter; ++n) {
					gradX += m_imageGauss[i + m][j + n] * filter[m + m_filterCenter][n + m_filterCenter];
					gradY += m_imageGauss[i + m][j + n] * filter[n + m_filterCenter][m + m_filterCenter];
				}
			}
			m_imageSobelY[i][j] = abs(gradY) > 65535 ? 65535 : abs(gradY);
			m_imageSobelX[i][j] = abs(gradX) > 65535 ? 65535 : abs(gradX);
			if (gradX == 0)
			{
				gradX = 0.00000000000000001;
			}
			m_pointDrection[i][j] = (atan(gradY / gradX) * 180) / pi;
			if(m_pointDrection[i][j]<0)
				m_pointDrection[i][j] += 180;
		}
	}
	SobelAmplitude();
}

//void MyCanny::SobelAmplitude()
//{
//	int height = m_imageSobelX.size();
//	int width = m_imageSobelX[0].size();
//	for (int i = m_filterCenter; i<height- m_filterCenter; i++)
//	{
//		for (int j = m_filterCenter; j<width- m_filterCenter; j++)
//		{
//			double ampXY = sqrt(m_imageSobelX[i][j] * m_imageSobelX[i][j] + m_imageSobelY[i][j] * m_imageSobelY[i][j]);
//			m_SobelAmpXY[i][j] = ampXY > 65535 ? 65535 : ampXY;
//		}
//	}
//}


//******************�Ǽ���ֵ����*************************
//��һ������imageInput�����Sobel�ݶ�ͼ��
//�ڶ�������imageOutPut������ľֲ�����ֵ����ͼ��
//���������pointDrection��ͼ����ÿ������ݶȷ�������ָ��
//*************************************************************
void MyCanny::LocalMaxValue(){
	if (m_SobelAmpXY.empty()) {
		cout << "MyCanny m_SobelAmpXY empty ��" << endl;
		return;
	}

	float g1 = 0, g2 = 0, g3 = 0, g4 = 0;                            //���ڽ��в�ֵ���õ������ص�����ֵ   
	double dTmp1 = 0.0, dTmp2 = 0.0;                           //�������������ص��ֵ�õ��ĻҶ����� 
	double dWeight = 0.0;                                    //��ֵ��Ȩ��  
	m_nonMaximum.assign(m_SobelAmpXY.begin(), m_SobelAmpXY.end());

	//�趨�ĽǶ���ֵΪ45��
	for (int i = 0; i < imageHeight ; i++)
	{
		for (int j = 0; j < imageWidth; j++)
		{
			if (m_SobelAmpXY[i][j] < 1e-6) {
				continue;
			}
			else {
				//178.5~180
				if (m_pointDrection[i][j] >= 178.5 && m_pointDrection[i][j]<180)
				{
					g1 = m_SobelAmpXY[i - 1][j - 1];
					g2 = m_SobelAmpXY[i][j - 1];
					g3 = m_SobelAmpXY[i][j + 1];
					g4 = m_SobelAmpXY[i + 1][j + 1];
					dWeight = m_imageSobelY[i][j] / m_imageSobelX[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
				//0~1.5
				else if (m_pointDrection[i][j] >= 0 && m_pointDrection[i][j]<1.5)
				{
					g1 = m_SobelAmpXY[i - 1][j + 1];
					g2 = m_SobelAmpXY[i][j + 1];
					g3 = m_SobelAmpXY[i][j - 1];
					g4 = m_SobelAmpXY[i + 1][j - 1];
					dWeight = m_imageSobelY[i][j] / m_imageSobelX[i][j];
					dTmp1 = g1*dWeight + (1 - dWeight)*g2;
					dTmp2 = g4*dWeight + (1 - dWeight)*g3;
				}
				else {
					m_nonMaximum[i][j] = 0;
					continue;
				}
			}
			if (m_nonMaximum[i][j] <= dTmp1 || m_nonMaximum[i][j] <= dTmp2)
				m_nonMaximum[i][j] = 0;
		}
	}
}

////˫��ֵ��ѡȡ
//void MyCanny::SetThreshold() {
//	int nHist[65536]{ 0 };//ֱ��ͼ
//	int nEdgeNum = 0;//���б�Ե�����Ŀ
//	int nMaxMag = 0;//����ݶȵķ�ֵ
//
//	int height = m_nonMaximum.size();
//	int width = m_nonMaximum[0].size();
//
//	for (int i = m_filterCenter; i < height- m_filterCenter; ++i) {
//		for (int j = m_filterCenter; j < width- m_filterCenter; ++j) {
//			int nindex = m_nonMaximum[i][j];
//			if (nindex > 0) {
//				nHist[nindex]++;
//				nEdgeNum++;
//				nMaxMag = max(nindex, nMaxMag);
//			}	
//		}
//	}
//
//	//����������ֵ ע�����ݶȵ���ֵ
//	double dRateHigh = 0.5;
//	double dRateLow = 0.5;
//	int nHightcount = (int)(dRateHigh*nEdgeNum + 0.5);
//	int count = 0;
//	nEdgeNum = nHist[1];
//	while ((nEdgeNum <= nHightcount) && (count < nMaxMag - 1)) {
//		count++;
//		nEdgeNum += nHist[count];
//	}
//	m_highThreshold = double(count);
//	m_lowThreshold = m_highThreshold*dRateLow;
//}


void MyCanny::DoubleThreshold()
{
	if (m_nonMaximum.empty()) {
		cout << "m_nonMaximum empty ��" << endl;
		return;
	}

	m_imageDouble.assign(m_nonMaximum.begin(), m_nonMaximum.end());
	m_imageEdgeStart.assign(m_nonMaximum.begin(), m_nonMaximum.end());

	//��ø���ֵ�͵���ֵ  ����ֵΪ50%������ֵΪ50%����ֵ
	SetThreshold();

	for (int i = 0; i < imageHeight; i++)
	{
		for (int j = 0; j < imageWidth; j++) 
		{
			m_imageEdgeStart[i][j] = 0;
			if (m_nonMaximum[i][j]>m_highThreshold)
			{
				m_imageDouble[i][j] = 65535;
				m_imageEdgeStart[i][j] = 65535;

			}
			else if (m_nonMaximum[i][j]<m_lowThreshold)
			{
				m_imageDouble[i][j] = 0;
			}
		}
	}
}

//void MyCanny::DoubleThresholdLink()
//{
//	m_imageLink.assign(m_imageDouble.begin(), m_imageDouble.end());
//	vector<pair<int,int>> flag;//������¼����˫��ֵ֮��ĵ�
//	int change = 0;//��¼�仯�ĵ�ĸ���
//	for (int i = m_filterCenter; i < m_imageLink.size() - m_filterCenter; i++)
//	{
//		for (int j = m_filterCenter; j < m_imageLink[0].size() - m_filterCenter; j++)
//		{
//			if (m_imageLink[i][j] > m_lowThreshold&&m_imageLink[i][j] <= m_highThreshold)
//			{
//				if (m_imageLink[i - 1][j - 1] == 65535 || m_imageLink[i - 1][j] == 65535 || m_imageLink[i - 1][j + 1] == 65535 ||
//					m_imageLink[i][j - 1] == 65535 || m_imageLink[i][j + 1] == 65535 || m_imageLink[i + 1][j - 1] == 65535 ||
//					m_imageLink[i + 1][j] == 65535 || m_imageLink[i + 1][j + 1] == 65535)
//				{
//					m_imageLink[i][j] = 65535;
//					change++;
//				}
//				else
//				{
//					//imageInput[i][j] = 0;
//					flag.push_back(make_pair(i, j));
//				}
//			}
//		}
//	}
//
//	while (change > 0&&!flag.empty()) {
//		change = 0;
//		vector<pair<int, int>> tempflag;
//		for (auto it:flag)
//		{
//			int i = it.first;
//			int j = it.second;
//			if (m_imageLink[i - 1][j - 1] == 65535 || m_imageLink[i - 1][j] == 65535 || m_imageLink[i - 1][j + 1] == 65535 ||
//				m_imageLink[i][j - 1] == 65535 || m_imageLink[i][j + 1] == 65535 || m_imageLink[i + 1][j - 1] == 65535 ||
//				m_imageLink[i + 1][j] == 65535 || m_imageLink[i + 1][j + 1] == 65535)
//			{
//				m_imageLink[i][j] = 65535;
//				change++;
//			}
//			else {
//				tempflag.push_back(it);
//			}
//		}
//		flag = tempflag;
//	}
//	for (auto it : flag) {
//		int i = it.first;
//		int j = it.second;
//		m_imageLink[i][j] = 0;
//	}
//}

void MyCanny::LengthThresholdLink() {
	if (m_nonMaximum.empty()) {
		cout << "MyCanny m_nonMaximum empty ��" << endl;
		return;
	}
	vector<vector<int>> flag0(imageHeight, vector<int>(imageWidth, 0));
	vector<vector<int>> flag1(flag0);
	m_imageLink.resize(imageHeight, vector<float>(imageWidth, 0));

	int T0 = 2;			//���ӳ�����ֵ
	int T1 = T0 * 2;

	vector<vector<float>> temp(m_nonMaximum);
	for (int j = 0; j < imageWidth; ++j) {
		int preEnd = 0;
		for (int i = 0; i < imageHeight; ++i) {
			if (flag0[i][j] == 1) continue;
			if (m_imageEdgeStart[i][j] == 65535) {
				int pre = preEnd;
				LinkHelp(flag0, T0, i, j,preEnd, temp);

				if (preEnd != pre && (i - pre <= T0 * 2)) {
					for (int m = pre + 1; m < i; ++m) {
						temp[m][j] = 65535;
					}
				}
			}
		}
	}

	//temp = m_imageLink;
	//m_imageLink.assign(height, vector<float>(width, 0));
	//for (int j = 0; j < width; ++j) {
	//	int preEnd = 0;
	//	for (int i = 0; i < height; ++i) {
	//		if (flag1[i][j] == 1) continue;
	//		if (m_imageLink[i][j] == 65535) {
	//			int pre = preEnd;
	//			LinkHelp(flag1, T1, i, j, preEnd, temp);

	//			if (preEnd != pre && (i - pre <= T1 * 2)) {
	//				for (int m = pre + 1; m < i; ++m) {
	//					m_imageLink[m][j] = 65535;
	//				}
	//			}
	//		}
	//	}
	//}
}

void MyCanny::LinkHelp(vector<vector<int>>& flag,int T,int i,int j,int& End,vector<vector<float>>& pImageIn) {
	int begin = i;
	int end = i;
	while (i < imageHeight && pImageIn[i][j] != 0) {
		flag[i][j] = 1;
		end = i;
		i++;
	}
	if (end - begin >= T) {
		End = end;
		for (int m = begin; m <= end; ++m) {
			m_imageLink[m][j] = 65535;
		}
	}
}

void MyCanny::CalEdgeStartToEnd() {
	for (int i = 0; i < imageWidth; ++i) {
		int start = 0, end = imageHeight - 1;
		while (start < end && m_imageEdgeStart[start][i] == 0) {
			start++;
		}
		while (end > start && m_imageEdgeStart[end][i] == 0) {
			end--;
		}
		m_startToEnd.push_back(make_pair(start, end));
	}
}

void MyCanny::GetArtifactEdge() {
	if (m_imageEdgeStart.empty()) {
		cout << "MyCanny m_imageEdgeStart empty ��" << endl;
		return;
	}
	CalEdgeStartToEnd();

	m_imageArtifact.resize(imageHeight, vector<float>(imageWidth, 0));
	float ratio = 0.07;		
	int start = 0;
	int end = 0;

	for (int i = 0; i < imageWidth; ++i) {
		start = m_startToEnd[i].first;
		end = m_startToEnd[i].second;
		if (end <= start) continue;
		int linknum = 0;
		//ͳ��ǿ��Ե��ʼ�����ֹ��֮����ڵ����ظ���
		for (int j = start; j <= end; ++j) {
			if (m_imageLink[j][i] == 65535)
				linknum++;
		}
		//�����ڵ����ظ�������һ������  �������ΪαӰ
		if (linknum >= (end - start)*ratio) {
			for (int j = start; j <= end; ++j) {
				m_imageArtifact[j][i] = 65535;
			}
		}
	}
}


void MyCanny::ArtifactCorrect(){
	if (m_srcImage.empty()) {
		cout << "MyCanny m_srcImage empty ��" << endl;
		return;
	}
	m_desImage.assign(m_srcImage.begin(), m_srcImage.end());

	int divide = 30;	//��ͼ������ֱ�����Ϊdivide�ν������
	int num = imageHeight / divide;
	vector<CPosition> bSlineNode;


	for (int k = 0; k < divide; ++k) {
		vector<pair<int,double>> artifact;	//αӰ���ĺ�����͸öε��о�ֵ

		for (int j = 0; j < imageWidth; ++j) {
			bool hasArtifact = false;	//αӰ�α�ʶ
			double	meanVal = 0;		//�ö� ÿ�е��о�ֵ

			for (int i = num*k; i < num*(k + 1); i++) {
				//�ж��ǲ���αӰ��
				if (m_imageArtifact[i][j] == 65535) hasArtifact = true;
				meanVal = meanVal + m_srcImage[i][j] / num;
			}

			if (!hasArtifact)	//����αӰ�Σ���ӵ�B�������������
				bSlineNode.push_back(CPosition(j, meanVal));
			else
				artifact.push_back(make_pair(j, meanVal));
		}


		ThreeOrderBSplineInterpolatePt(bSlineNode);

		//У��αӰ��������ֵ
		for (int m = 0; m < artifact.size(); ++m) {
			double gainFactor = 0;	//αӰ����������
			int index_x = artifact[m].first;
			gainFactor = bSlineNode[index_x].y / artifact[m].second;

			for (int n = num*k; n < num*(k + 1); n++) {
				m_desImage[n][index_x] = gainFactor*m_srcImage[n][index_x];
			}	
		}
	}

}


//B�������

//================================================================
// �������ܣ� ����B�������,�ڽڵ�֮����Ȳ���ָ��������
// ��������� *pt �����������У�ִ����ɺ󣬻ᱻ�滻���µ����ݵ�
//            Num���ڵ�����
//            *InsertNum: �ڵ�֮����Ҫ����ĵ����ָ�� 
// ����ֵ��   �޷���ֵ
//
//=================================================================
void MyCanny::ThreeOrderBSplineInterpolatePt(vector<CPosition>& pt){
	if (pt.size() < 5) {
		return;
	}

	int Num = pt.size();
	vector<CPosition> temp;
	double x0 = 2 * pt[0].x - pt[1].x;
	double y0= 2 * pt[0].y - pt[1].y;
	temp.push_back(CPosition(x0, y0));		  //  �������ӳ����ϵ������Ϊ�׵�
	for (int i = 0; i<Num; i++)
		temp.push_back(pt[i]);

	x0 = 2 * temp[Num].x - temp[Num - 1].x;
	y0 = 2 * temp[Num].y - temp[Num - 1].y;
	temp.push_back(CPosition(x0, y0));		  //  �������ӳ����ϵ������Ϊβ��

	CPosition NodePt1, NodePt2, NodePt3, NodePt4;
	double t;

	pt.clear();

	int index = 0;
	for (int i = 0; i<Num - 1; i++)                          //  ÿ���߶ξ��Ȳ����
	{
		NodePt1 = temp[i]; NodePt2 = temp[i + 1]; NodePt3 = temp[i + 2]; NodePt4 = temp[i + 3];
		int InsertNum = temp[i + 1].x - temp[i].x;
		double dt = 1.0 / InsertNum;

		for (int j = 0; j<InsertNum ; j++)
		{
			t = dt*j;
			x0 = index;
			y0 = F03(t)*NodePt1.y + F13(t)*NodePt2.y + F23(t)*NodePt3.y + F33(t)*NodePt4.y;
			pt.push_back(CPosition(x0, y0));
			index++;
		}

		if (i == Num - 2) {              //  ���һ��β��
			t = 1;
			x0 = pt.back().x + 1;
			y0 = F03(t)*NodePt1.y + F13(t)*NodePt2.y + F23(t)*NodePt3.y + F33(t)*NodePt4.y;
			pt.push_back(CPosition(x0, y0));
		}
	}
}

//================================================================
// �������ܣ� ��������������
//
//================================================================
double MyCanny::F03(double t)
{
	return 1.0 / 6 * (-t*t*t + 3 * t*t - 3 * t + 1);
}
double MyCanny::F13(double t)
{
	return 1.0 / 6 * (3 * t*t*t - 6 * t*t + 4);
}
double MyCanny::F23(double t)
{
	return 1.0 / 6 * (-3 * t*t*t + 3 * t*t + 3 * t + 1);
}
double MyCanny::F33(double t)
{
	return 1.0 / 6 * t*t*t;
}
