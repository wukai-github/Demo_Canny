#include <vector>
using namespace std;

//B������
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

	//���ò���
	inline void setImageWidth(int width) { imageWidth = width; };
	inline void setImageHight(int height) { imageHeight = height; };
	inline void setGuassKernelSize(int size) { guassKernelSize = size; };
	void setSrcImage(unsigned short* srcImage);
	//��ȡ����
	inline const int getImageWidth() const { return imageWidth; };
	inline const int getImageHeight() const { return imageHeight; };

	//��ȡͼ��
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
	void Gaussianconv();//��˹�˲�
	virtual void SobelGradDirction();	//����ͼ���ݶȣ������ݶȼ���ͼ���Ե��ֵ��Ƕ�
	virtual void LocalMaxValue();		//�Ǽ���ֵ���ƣ���Եϸ����
	virtual void DoubleThreshold();		//˫��ֵ����
	void DoubleThresholdLink();	//�м����ش�����Ե����

	void SobelAmplitude();		//����ͼ���Ե��ֵ��X,Y����ĵ��ӣ�
	void GetGaussianKernel(vector<vector<float>>& Kernel, const float sigmma = 1);		//��ȡ��˹��
	void SetThreshold();		//�趨˫��ֵ

protected:
	int imageWidth;
	int imageHeight;
	int guassKernelSize;

	vector<vector<float>> m_srcImage;		//ԭʼͼ��
	vector<vector<float>> m_imageGauss;		//��˹ģ��ͼ��
	vector<vector<float>> m_imageSobelX;	//X�����ݶ�ֵ
	vector<vector<float>> m_imageSobelY;	//Y�����ݶ�ֵ
	vector<vector<float>> m_pointDrection;	//�Ƕ�
	vector<vector<float>> m_SobelAmpXY;		//ͼ���ݶ�ֵ
	vector<vector<float>> m_nonMaximum;		//�Ǽ���ֵ����ͼ��
	vector<vector<float>> m_imageDouble;	//˫��ֵ������ͼ��
	vector<vector<float>> m_imageLink;		//��Ե���Ӻ��ͼ��

	double m_lowThreshold = 0;	//˫��ֵ�ĵ���ֵ
	double m_highThreshold = 0;	//˫��ֵ�ĸ���ֵ

};


class MyCanny :public CannyBase{
	/**************************��״αӰʶ��***********************/

public:
	MyCanny(unsigned short* srcImage, int width, int height, int size);
public:
	//void Gaussianconv(unsigned short* pImageIn, int width, int height, const int size = 5);//��˹�˲�
	virtual void SobelGradDirction() override;	//����ͼ���ݶȣ������ݶȼ���ͼ���Ե��ֵ��Ƕ�
	virtual void LocalMaxValue() override;		//�Ǽ���ֵ���ƣ���Եϸ����
	virtual void DoubleThreshold() override;		//˫��ֵ����
	//void DoubleThresholdLink();	//�м����ش�����Ե����
	void LengthThresholdLink();	//�趨������ֵ��������
	void GetArtifactEdge();		//��ȡ��״αӰ

	//vector<vector<float>> m_srcImage;		//ԭʼͼ��
	//vector<vector<float>> m_imageGauss;		//��˹ģ��ͼ��
	//vector<vector<float>> m_imageSobelX;	//X�����ݶ�ֵ
	//vector<vector<float>> m_imageSobelY;	//Y�����ݶ�ֵ
	//vector<vector<float>> m_pointDrection;	//�Ƕ�
	//vector<vector<float>> m_SobelAmpXY;		//ͼ���ݶ�ֵ
	//vector<vector<float>> m_nonMaximum;		//�Ǽ���ֵ����ͼ��
	//vector<vector<float>> m_imageDouble;	//˫��ֵ������ͼ��
	vector<vector<float>> m_imageEdgeStart;	//��Ե��ʼ��
	//vector<vector<float>> m_imageLink;		//��Ե���Ӻ��ͼ��
	vector<vector<float>> m_imageArtifact;	//��״αӰ


private:
	//void GetGaussianKernel(vector<vector<float>>& Kernel,const int KernelSize = 5, const float sigmma = 1);//��ȡ��˹��
	//void SobelAmplitude();		//����ͼ���Ե��ֵ��X,Y����ĵ��ӣ�
	//void SetThreshold();			//�趨˫��ֵ
	void CalEdgeStartToEnd();		//��ÿ��ǿ��Ե�����Ϻ����µ��������������

	void LinkHelp(vector<vector<int>>& flag, int T,int i,int j,int& End, vector<vector<float>>& pImageIn);			//����

	int m_filtersize = 7;		//�˲�ģ��ߴ�
	int m_filterCenter = m_filtersize / 2;
	
	//double m_lowThreshold = 0;	//˫��ֵ�ĵ���ֵ
	//double m_highThreshold = 0;	//˫��ֵ�ĸ���ֵ

	vector<pair<int, int>> m_startToEnd;	//��¼ÿ��ǿ��Ե�����Ϻ����µ�������������� 


	/*************************��״αӰУ��**************************************/
public:
	void ArtifactCorrect();		//У����״αӰ
	vector<vector<float>> m_desImage;		//���ͼ��

/**************************����B�������**********************************/
public:
	//void ThreeOrderBSplineSmooth(CPosition *pt, int Num);		B����ƽ��
	void ThreeOrderBSplineInterpolatePt(vector<CPosition>& pt);	//B������� ��ֵ

private:
	double F03(double t);
	double F13(double t);
	double F23(double t);
	double F33(double t);

};
