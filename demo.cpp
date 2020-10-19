#include <iostream>
#include "BinFile.h"
#include "opencv2/core/core.hpp"
#include "opencv2/opencv.hpp"  
#include "Canny.h"

using namespace cv;

MyCanny m_canny;

void Canny(unsigned short *imageData,int imageHeight,int imageWidth) {
	//ԭʼͼ��
	Mat srcImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			srcImage.at<unsigned short>(i, j) = imageData[i*imageWidth + j];
		}
	}
	namedWindow("srcImage", 0);
	cvResizeWindow("srcImage", imageWidth / 3, imageHeight / 3);
	imshow("srcImage", srcImage);

	//��˹���
	//m_canny.Gaussianconv();
	Mat gauImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			gauImage.at<unsigned short>(i, j) = m_canny.getGuassImage()[i][j];
		}
	}
	namedWindow("gauImage", 0);
	cvResizeWindow("gauImage", imageWidth / 3, imageHeight / 3);
	imshow("gauImage", gauImage);

	//ͼ���ݶ�
	m_canny.SobelGradDirction();
	Mat sobelImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			sobelImage.at<unsigned short>(i, j) = m_canny.getSobelImage()[i][j];
		}
	}
	namedWindow("sobelImage", 0);
	cvResizeWindow("sobelImage", imageWidth / 3, imageHeight / 3);
	imshow("sobelImage", sobelImage);

	//�Ǽ���ֵ����
	m_canny.LocalMaxValue();
	Mat nonMaxImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			nonMaxImage.at<unsigned short>(i, j) = m_canny.getNonMaxImage()[i][j];
		}
	}
	namedWindow("nonMaxImage", 0);
	cvResizeWindow("nonMaxImage", imageWidth / 3, imageHeight / 3);
	imshow("nonMaxImage", nonMaxImage);

	//˫��ֵ����
	m_canny.DoubleThreshold();
	//���ӱ�Ե
	m_canny.DoubleThresholdLink();
	//m_canny.LengthThresholdLink();
	//�õ�αӰͼ��
	m_canny.GetArtifactEdge();
	Mat linkImage(imageHeight, imageWidth, CV_16UC1);
	Mat ArtifactImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			linkImage.at<unsigned short>(i, j) = m_canny.getLinkImage()[i][j];
			ArtifactImage.at<unsigned short>(i, j) = m_canny.m_imageArtifact[i][j];
		}
	}
	namedWindow("linkImage", 0);
	cvResizeWindow("linkImage", imageWidth / 3, imageHeight / 3);
	imshow("linkImage", linkImage);
	namedWindow("ArtifactImage", 0);
	cvResizeWindow("ArtifactImage", imageWidth / 3, imageHeight / 3);
	imshow("ArtifactImage", ArtifactImage);

	waitKey(0);
}

int main()
{
    char *pathname = "D:\\work\\ZZ\\2020-01-02-15-54-16\\Filter Data\\filt_0000.bin";
	//char *pathname = "D:\\work\\ZZ\\2020-01-02-15-54-16\\Projection Data\\proj_0000.bin";
	
    unsigned short *imageData = nullptr;
    int imageWidth=0;
    int imageHeight=0;
    double dmax=0;
    double dmin=0;

	//��ȡBinͼ��
    bool isRead = ReadBin(pathname, imageData, imageWidth, imageHeight, dmax, dmin);
	if (!isRead) {
		cout << "ͼ�񲻴��ڣ�" << endl;
		return 0;
	}
	
	//Canny(imageData, imageHeight, imageWidth);
	//ԭʼͼ��
	Mat srcImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			srcImage.at<unsigned short>(i, j) = imageData[i*imageWidth + j];
		}
	}
	namedWindow("srcImage", 0);
	cvResizeWindow("srcImage", imageWidth / 3, imageHeight / 3);
	imshow("srcImage", srcImage);

	//��˹���
	//m_canny.Gaussianconv();
	Mat gauImage(imageHeight, imageWidth, CV_16UC1);
	for (int i = 0; i<imageHeight; ++i) {
		for (int j = 0; j<imageWidth; ++j) {
			gauImage.at<unsigned short>(i, j) = m_canny.getGuassImage()[i][j];
		}
	}
	namedWindow("gauImage", 0);
	cvResizeWindow("gauImage", imageWidth / 3, imageHeight / 3);
	imshow("gauImage", gauImage);

	//���з������ͼ��
	vector<int> integral(imageWidth, 0);
	Mat integralImage(imageHeight, imageWidth, CV_16UC1);
	long sum = 0;	//ÿ�лҶ�ֵ���ܺͣ�ӳ�䵽0-imageheight
	for (int i = 0; i < imageWidth; ++i) {
		sum = 0;
		for (int j = 0; j < imageHeight; ++j) {
			//sum += gauImage.at<unsigned short>(j, i);
			sum += srcImage.at<unsigned short>(j, i);
		}
		sum = sum / 65535;
		integral[i] = sum;
		for (int k = imageHeight-1; k >= 0; --k) {
			if(k>= imageHeight-sum)
				integralImage.at<unsigned short>(k, i) = 0;
			else
				integralImage.at<unsigned short>(k, i) = 65535;
		}
	}
	namedWindow("integralImage", 0);
	cvResizeWindow("integralImage", imageWidth / 3, imageHeight / 3);
	imshow("integralImage", integralImage);
	
	//һ��΢�ִ�����ͼ��
	vector<int> diffData(imageWidth, 0);		//��¼ÿ�е�һ��΢��ֵ
	int maxdiff = 0;	//һ��΢���е����ֵ
	Mat diffImage(imageHeight, imageWidth, CV_16UC1);
	int diff = 0;		//һ��΢��ֵ
	for (int j = 0; j < imageWidth; ++j) {
		//�����һ��
		if (j == 0) {
			for (int i = 0; i < imageHeight; ++i) {
				diffImage.at<unsigned short>(i, j) = 65535;
			}
		}
		else {
			diff = abs(integral[j] - integral[j - 1]);
			diffData[j] = diff;
			maxdiff = max(maxdiff, diff);
			for (int k = imageHeight - 1; k >= 0; --k) {
				if (k >= imageHeight - diff)
					diffImage.at<unsigned short>(k, j) = 0;
				else
					diffImage.at<unsigned short>(k, j) = 65535;
			}
		}
	}
	namedWindow("diffImage", 0);
	cvResizeWindow("diffImage", imageWidth / 3, imageHeight / 3);
	imshow("diffImage", diffImage);

	//����һ��΢��ֵ��ֱ��ͼ�����ƽ��ֵ�;����� ���趨αӰ��ֵ
	int histWidth = maxdiff + 1;
	vector<int> hist(histWidth, 0);		//ֱ��ͼ����
	double average = 0;		//ƽ��ֵ
	double variance = 0;	//������
	int coefficient = 3.0;	//����ϵ��  threshold=average+coefficient*variance
	double threshold = 0;	
	Mat histImage(imageHeight, histWidth, CV_16UC1);	//ֱ��ͼ����
	//��ֱ��ͼ����ֵ
	for (int i = 0; i < imageWidth; ++i) {
		int index = diffData[i];
		average += diffData[i];
		hist[index]++;
	}
	for (int j = 0; j < histWidth; ++j) {
		for (int i = imageHeight - 1; i >= 0; --i) {
			if (i >= imageHeight - hist[j])
				histImage.at<unsigned short>(i, j) = 0;
			else
				histImage.at<unsigned short>(i, j) = 65535;
		}
	}
	average /= imageWidth;
	for (int i = 0; i < imageWidth; ++i) {
		variance += pow((diffData[i] - average), 2);
	}
	variance = sqrt(variance / (imageWidth - 1));
	threshold = average + coefficient*variance;
	//namedWindow("histImage", 0);
	//cvResizeWindow("histImage", histWidth / 3, imageHeight / 3);
	//imshow("histImage", histImage);

	/*һ��΢�������ֵ���������ֵ����ΪαӰ��*/
	//Ѱ��һ��΢������Ķ�Ӧ��ֵ
	//int threshIndex = maxdiff;		//һ��΢������Ķ�Ӧ��ֵ
	//for (int i = maxdiff; i >= 0; --i) {
	//	if (hist[i] > threshold){
	//		threshIndex = i;
	//		break;
	//	}
	//}
	//if (threshIndex == maxdiff) return 0;
	vector<bool> artifactIndex(imageWidth, false);		//αӰ��
	for (int i = 0; i < imageWidth; ++i) {
		if (diffData[i] > threshold)
			artifactIndex[i] = true;
	}


	//αӰͼ��
	Mat  artifactImage;
	srcImage.copyTo(artifactImage);
	for (int j = 0; j < imageWidth; ++j) {
		if(artifactIndex[j])
			for (int i = 0; i < imageHeight; ++i) {
				artifactImage.at<unsigned short>(i, j) = 0;
			}
	}
	namedWindow("artifactImage", 0);
	cvResizeWindow("artifactImage", imageWidth / 3, imageHeight / 3);
	imshow("artifactImage", artifactImage);

	//αӰУ��  �������Ʋ�ֵ
	Mat dstImage(imageHeight, imageWidth, CV_16UC1);	//У�����ͼ��
	srcImage.copyTo(dstImage);
	for (int j = 2; j < imageWidth - 2; ++j) {
		if (artifactIndex[j]) {
			int start = j;	//����αӰ�е���ʼ������
			while (j < imageWidth&&artifactIndex[j + 1])	++j;
			int end = j;	//����αӰ�е��յ�������
			bool nextFlag = false;		//�ж��ұ��Ƿ�������������
			if (!artifactIndex[j + 1] && !artifactIndex[j + 2])
				nextFlag = true;

			//�Ե���αӰ���ݵ��������ڴ���ͶӰֵ������������  У��
			//double d = 0;		//��ֵ��Ȩֵ
			//if (end == start) {
			//	for (int i = 0; i < imageHeight; ++i) {
			//	/*	if(start>3)
			//		dstImage.at<unsigned short>(i, start - 1) =
			//		2 * dstImage.at<unsigned short>(i, start - 2) - dstImage.at<unsigned short>(i, start - 3);
			//		if (nextFlag&&end < imageWidth - 3)
			//		dstImage.at<unsigned short>(i, end + 1) =
			//		2 * dstImage.at<unsigned short>(i, end + 2) - dstImage.at<unsigned short>(i, end + 3);*/
			//		for (int k = start; k <= end; ++k) {
			//			d = (double)(k - start + 1) / (end - start + 2);
			//			dstImage.at<unsigned short>(i, k) =
			//				(1 - d)*dstImage.at<unsigned short>(i, start - 1) + d*dstImage.at<unsigned short>(i, end + 1);
			//		}
			//	}
			//}

			//�Զ���αӰ���ݵ��ø�˹����������ֵ����  У��
			for (int i = 0; i < imageHeight; ++i) {
				for (int k = start; k <= end; ++k) {
					dstImage.at<unsigned short>(i, k) = gauImage.at<unsigned short>(i, k);
				}
			}
			/*else {
				for (int i = 0; i < imageHeight; ++i) {
					for (int k = start; k <= end; ++k) {
						dstImage.at<unsigned short>(i, k) = gauImage.at<unsigned short>(i, k);
					}
				}
			}*/
		}
	}
	namedWindow("dstImage", 0);
	cvResizeWindow("dstImage", imageWidth / 3, imageHeight / 3);
	imshow("dstImage", dstImage);

	float *CorrectImage = new float[imageWidth*imageHeight];
	double dVmin = 65535, dVmax = 0;
	for (int i = 0; i < imageHeight; ++i) {
		for (int j = 0; j < imageWidth; ++j) {
			CorrectImage[i*imageWidth + j] = (float)dstImage.at<unsigned short>(i, j);
			dVmin = min(dVmin, (double)CorrectImage[i*imageWidth + j]);
			dVmax = max(dVmax, (double)CorrectImage[i*imageWidth + j]);
		}
	}
	WriteBin("CorrectImage.bin", CorrectImage, imageWidth, imageHeight, dVmin, dVmax);

	//char *pathname1 = "CorrectImage.bin";

	//unsigned short *imageData1 = nullptr;
	//int imageWidth1 = 0;
	//int imageHeight1 = 0;
	//double dmax1 = 0;
	//double dmin1 = 0;

	////��ȡBinͼ��
	//bool isRead1 = ReadBin(pathname1, imageData1, imageWidth1, imageHeight1, dmax1, dmin1);

	////Canny(imageData, imageHeight, imageWidth);
	////ԭʼͼ��
	//Mat srcImage1(imageHeight, imageWidth, CV_16UC1);
	//for (int i = 0; i<imageHeight; ++i) {
	//	for (int j = 0; j<imageWidth; ++j) {
	//		srcImage1.at<unsigned short>(i, j) = imageData1[i*imageWidth + j];
	//	}
	//}
	//namedWindow("srcImage1", 0);
	//cvResizeWindow("srcImage1", imageWidth1 / 3, imageHeight1 / 3);
	//imshow("srcImage1", srcImage1);

	waitKey(0);
    return 0;
}
