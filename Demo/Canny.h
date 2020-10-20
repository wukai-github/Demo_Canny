#include <vector>
using namespace std;

//B样条点
typedef struct tagPosition
{
	int  x;
	double  y;
	tagPosition(int _x, double _y) { x = _x; y = _y; }
	tagPosition() {};
	bool operator==(const tagPosition & pt) { return (x == pt.x && y == pt.y); }
} CPosition;

class CannyBase {
public:
	CannyBase() {};
	CannyBase(unsigned short* srcImage, int width, int height, int size);
	~CannyBase() {};

	//设置参数
	inline void setImageWidth(int width) { imageWidth = width; };
	inline void setImageHight(int height) { imageHeight = height; };
	inline void setGuassKernelSize(int size) { guassKernelSize = size; };
	void setSrcImage(unsigned short* srcImage);
	//获取参数
	inline const int getImageWidth() const { return imageWidth; };
	inline const int getImageHeight() const { return imageHeight; };

	//获取图像
	inline const vector<vector<float>> getSrcImage() const { return m_srcImage; };
	inline const vector<vector<float>> getGuassImage() const { return m_imageGauss; };
	inline const vector<vector<float>> getSobelImage() const { return m_SobelAmpXY; };
	inline const vector<vector<float>> getSobelXImage() const { return m_imageSobelX; };
	inline const vector<vector<float>> getSobelYImage() const { return m_imageSobelY; };
	inline const vector<vector<float>> getNonMaxImage() const { return m_nonMaximum; };
	inline const vector<vector<float>> getDoubleImage() const { return m_imageDouble; };
	inline const vector<vector<float>> getLinkImage() const { return m_imageLink; };

	void CannyEdgeDetect();
	

protected:
	void Gaussianconv();//高斯滤波
	virtual void SobelGradDirction();	//计算图像梯度，根据梯度计算图像边缘幅值与角度
	virtual void LocalMaxValue();		//非极大值抑制（边缘细化）
	virtual void DoubleThreshold();		//双阈值处理
	void DoubleThresholdLink();	//中间像素处理及边缘链接

	void SobelAmplitude();		//计算图像边缘幅值（X,Y方向的叠加）
	void GetGaussianKernel(vector<vector<float>>& Kernel, const float sigmma = 1);		//获取高斯核
	void SetThreshold();		//设定双阈值

protected:
	int imageWidth;
	int imageHeight;
	int guassKernelSize;

	vector<vector<float>> m_srcImage;		//原始图像
	vector<vector<float>> m_imageGauss;		//高斯模糊图像
	vector<vector<float>> m_imageSobelX;	//X方向梯度值
	vector<vector<float>> m_imageSobelY;	//Y方向梯度值
	vector<vector<float>> m_pointDrection;	//角度
	vector<vector<float>> m_SobelAmpXY;		//图像梯度值
	vector<vector<float>> m_nonMaximum;		//非极大值抑制图像
	vector<vector<float>> m_imageDouble;	//双阈值处理后的图像
	vector<vector<float>> m_imageLink;		//边缘连接后的图像

	double m_lowThreshold = 0;	//双阈值的低阈值
	double m_highThreshold = 0;	//双阈值的高阈值

};


class MyCanny :public CannyBase{
	/**************************环状伪影识别***********************/

public:
	MyCanny(unsigned short* srcImage, int width, int height, int size);
public:
	//void Gaussianconv(unsigned short* pImageIn, int width, int height, const int size = 5);//高斯滤波
	virtual void SobelGradDirction() override;	//计算图像梯度，根据梯度计算图像边缘幅值与角度
	virtual void LocalMaxValue() override;		//非极大值抑制（边缘细化）
	virtual void DoubleThreshold() override;		//双阈值处理
	//void DoubleThresholdLink();	//中间像素处理及边缘链接
	void LengthThresholdLink();	//设定长度阈值进行连接
	void GetArtifactEdge();		//获取环状伪影

	//vector<vector<float>> m_srcImage;		//原始图像
	//vector<vector<float>> m_imageGauss;		//高斯模糊图像
	//vector<vector<float>> m_imageSobelX;	//X方向梯度值
	//vector<vector<float>> m_imageSobelY;	//Y方向梯度值
	//vector<vector<float>> m_pointDrection;	//角度
	//vector<vector<float>> m_SobelAmpXY;		//图像梯度值
	//vector<vector<float>> m_nonMaximum;		//非极大值抑制图像
	//vector<vector<float>> m_imageDouble;	//双阈值处理后的图像
	vector<vector<float>> m_imageEdgeStart;	//边缘起始点
	//vector<vector<float>> m_imageLink;		//边缘连接后的图像
	vector<vector<float>> m_imageArtifact;	//环状伪影


private:
	//void GetGaussianKernel(vector<vector<float>>& Kernel,const int KernelSize = 5, const float sigmma = 1);//获取高斯核
	//void SobelAmplitude();		//计算图像边缘幅值（X,Y方向的叠加）
	//void SetThreshold();			//设定双阈值
	void CalEdgeStartToEnd();		//求每列强边缘的最上和最下的两个点的纵坐标

	void LinkHelp(vector<vector<int>>& flag, int T,int i,int j,int& End, vector<vector<float>>& pImageIn);			//连接

	int m_filtersize = 7;		//滤波模板尺寸
	int m_filterCenter = m_filtersize / 2;
	
	//double m_lowThreshold = 0;	//双阈值的低阈值
	//double m_highThreshold = 0;	//双阈值的高阈值

	vector<pair<int, int>> m_startToEnd;	//记录每列强边缘的最上和最下的两个点的纵坐标 


	/*************************环状伪影校正**************************************/
public:
	void ArtifactCorrect();		//校正环状伪影
	vector<vector<float>> m_desImage;		//输出图像

/**************************三次B样条拟合**********************************/
public:
	//void ThreeOrderBSplineSmooth(CPosition *pt, int Num);		B样条平滑
	void ThreeOrderBSplineInterpolatePt(vector<CPosition>& pt);	//B样条拟合 插值

private:
	double F03(double t);
	double F13(double t);
	double F23(double t);
	double F33(double t);

};
