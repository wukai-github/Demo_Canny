//
// Created by zhangxm on 2020/4/2.
//

#ifndef BINFILE_H
#define BINFILE_H

#define PIXEL_MAX 65535

typedef struct BIN_INFOHEAD
{
    // *.BIN file header struct
    char s[492]; // Reserved
    float min; // Minimal value of data
    float max; // Maximal value of data
    int width; // Width of data
    int height; // Height of data
    int depth; // Depth of data (slices)
} BIN_HEADER;

bool ReadBin(char* PathName,unsigned short*& ImageData,int& width,int& height, double& dmax, double& dmin);
bool WriteBin(char* name, float *p, int width, int height, double dVmin, double dVmax);

#endif //BINFILE_H
