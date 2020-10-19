//
// Created by zhangxm on 2020/4/2.
//
#include "BinFile.h"
#include <iostream>

bool ReadBin(char* PathName, unsigned short*& ImageData, int& width, int& height, double& dmax, double& dmin) {
    //读模式打开图像文件
    FILE* fp = nullptr;
    if ((fp = fopen(PathName, "rb")) == nullptr)
        return false;

    int HdrSize = sizeof(BIN_HEADER);
    fseek(fp, -HdrSize, SEEK_END);

    BIN_HEADER bin_hdr;
    if (fread(&bin_hdr, sizeof(BIN_HEADER), 1, fp) != 1) {
        fclose(fp);
        return false;
    }

    //读取BIN文件信息头数据
    width = bin_hdr.width;
    height = bin_hdr.height;
    dmax = bin_hdr.max;
    dmin = bin_hdr.min;
    double dSacle = (dmax - dmin) / PIXEL_MAX;

    fseek(fp, 0, SEEK_SET);
    ImageData = new unsigned short[width * height];
    float* pTempImage = new float[width * height];
    if (pTempImage == nullptr || ImageData == nullptr) {
        fclose(fp);
        return false;
    }

    fread(pTempImage, width * height * sizeof(float), 1, fp);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (dmax > 2000)
                ImageData[i * width + j] =
                        static_cast<unsigned short>(pTempImage[i * width + j]);
            else
                ImageData[i * width + j] =
                        static_cast<unsigned short>(((pTempImage[i * width + j] - dmin) / dSacle));
        }
    }

    if (pTempImage != nullptr) {
        delete[] pTempImage;
        pTempImage = nullptr;
    }

    fclose(fp);
    return true;
}

bool WriteBin(char* name, float *p, int width, int height, double dVmin, double dVmax)
{
	FILE* f = fopen(name, "w+b");	
	if (f == NULL) return false;
	float *pf, *buf = NULL;
	bool bResult = false;
	int size = width * sizeof(float);
	BIN_HEADER hdr;
	if (NULL == (buf = (float*)malloc(width * sizeof(float))))
	{
		printf("IDS_NO_MEMORY");
		goto out;
	}
	memset(hdr.s, ' ', sizeof(hdr.s));
	hdr.width = width;
	hdr.height = height;
	hdr.depth = 1;
	hdr.min = dVmin;
	hdr.max = dVmax;


	for (int j = 0; j<height; j++)
	{
		pf = buf;
		for (int i = 0; i<width; i++)
			*pf++ = *p++;
		if (1 != fwrite(buf, size, 1, f))
			goto out;
	}
	if (1 != fwrite(&hdr, sizeof(hdr), 1, f))
		goto out;
	bResult = true;
out:
	if (buf)
		free(buf);
	fclose(f);

	return bResult;
}