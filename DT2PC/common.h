#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <assert.h>
#include "Math3.h"
#include "json.h"

enum YUVformat
{
	YUV444,
	YUV420
};

class Parameter
{
public:
	std::string outputPlyFolder;
	std::string jsonFilePath;

	size_t frameNumber;
	std::string ContentName;

	//view ����
	std::vector<std::string> viewNames;
	std::vector<uint16_t> widths;
	std::vector<uint16_t> heights;
	std::vector<Vector3<double>> Ts; //�����ӵ�ƽ����������json�ļ������õ�������������
	std::vector<Vector3<double>> Rs; //�����ӵ���ת������json�ļ��ļ������õ�������������
	std::vector<double> Rnears;
	std::vector<double> Rfars;
	std::vector<std::string> projections;

	//����ͼ��Ϣ
	std::vector<std::string> textureFileNames;
	std::vector<uint16_t> texturebitDepths;
	std::vector<YUVformat> textureYUVFormats;
	//���ͼ��Ϣ
	std::vector<std::string> depthFileNames;
	std::vector<YUVformat> depthYUVFormats;
	std::vector<uint16_t> depthbitDepths;

	//ȫ��ͼ��Ϣ
	std::vector<std::pair<double, double>> horRanges;
	std::vector<std::pair<double, double>> verRanges;
	//͸��ͼ��Ϣ
	std::vector<double> deltaXs;
	std::vector<double> deltaYs;
	std::vector<double> fxs;
	std::vector<double> fys;

	void Parse(int argc, char* argv[])
	{
		for (int iArgIdx = 1; iArgIdx < argc; iArgIdx += 2)
		{
			if (strcmp("-h", argv[iArgIdx]) == 0) { PrintHelp(); exit(1); }
			else if (strcmp("-o", argv[iArgIdx]) == 0) outputPlyFolder = argv[iArgIdx + 1];
			else if (strcmp("-j", argv[iArgIdx]) == 0) jsonFilePath = argv[iArgIdx + 1];
			else { printf("���棺����ĵ����������%s!\n", argv[iArgIdx]);  PrintHelp();}
		}
		//��������в���
		if (outputPlyFolder.empty()) { printf("���棺δ��������ļ�����·��!\n"); outputPlyFolder = "."; }
		if (jsonFilePath.empty()) { printf("����δ����json�ļ�·��!\n"); PrintHelp(); exit(0); }

		//����json�ļ�
		Json::Value value;
		Json::Reader jsonReader;
		std::ifstream ifp;
		ifp.open(jsonFilePath);
		if (!ifp.is_open())
		{
			printf("����json�ļ���ʧ�ܣ�\n");
			exit(0);
		}

		jsonReader.parse(ifp, value);
		if (value["Content_name"].isNull())
		{
			printf("����:json�ļ���ȱ��\"Content_name\"���ݣ�\n");
			exit(0);
		}
		ContentName = value["Content_name"].asString();

		if (value["Frames_number"].isNull())
		{
			printf("����:json�ļ���ȱ��\"Frames_number\"���ݣ�\n");
			exit(0);
		}
		frameNumber = value["Frames_number"].asUInt64();

		if (value["cameras"].isNull())
		{
			printf("����:json�ļ���ȱ��\"cameras\"���ݣ�\n");
			exit(0);
		}
		for (int i = 0; i < value["cameras"].size(); i++)
		{
			if (value["cameras"][i]["Name"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Name\"���ݣ�\n", i + 1);
				exit(0);
			}
			viewNames.push_back(value["cameras"][i]["Name"].asString());

			if (value["cameras"][i]["textureFilePath"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"textureFilePath\"���ݣ�\n", i + 1);
				exit(0);
			}
			textureFileNames.push_back(value["cameras"][i]["textureFilePath"].asString());

			if (value["cameras"][i]["depthFilePath"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"depthFilePath\"���ݣ�\n", i + 1);
				exit(0);
			}
			depthFileNames.push_back(value["cameras"][i]["depthFilePath"].asString());

			if (value["cameras"][i]["Position"].isNull() || (value["cameras"][i]["Position"].size() != 3))
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Position\"���ݻ���\"Position\"���ݴ���\n", i + 1);
				exit(0);
			}
			Vector3<double> T;
			for (int j = 0; j < value["cameras"][i]["Position"].size(); j++)
			{
				T[j] = value["cameras"][i]["Position"][j].asDouble();
			}
			Ts.push_back(T);

			if (value["cameras"][i]["Rotation"].isNull() || (value["cameras"][i]["Rotation"].size() != 3))
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Rotation\"���ݻ���\"Rotation\"���ݴ���\n", i + 1);
				exit(0);
			}
			Vector3<double> R;
			for (int j = 0; j < value["cameras"][i]["Rotation"].size(); j++)
			{
				R[j] = value["cameras"][i]["Rotation"][j].asDouble();
			}
			Rs.push_back(R);

			if (value["cameras"][i]["Resolution"].isNull() || (value["cameras"][i]["Resolution"].size() != 2))
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Resolution\"���ݻ���\"Resolution\"���ݴ���\n", i + 1);
				exit(0);
			}
			uint16_t Res[2];
			for (int j = 0; j < value["cameras"][i]["Resolution"].size(); j++)
			{
				Res[j] = value["cameras"][i]["Resolution"][j].asUInt();
			}
			widths.push_back(Res[0]);
			heights.push_back(Res[1]);

			if (value["cameras"][i]["Projection"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Projection\"���ݣ�\n", i + 1);
				exit(0);
			}
			projections.push_back(value["cameras"][i]["Projection"].asString());

			if (value["cameras"][i]["Depth_range"].isNull() || (value["cameras"][i]["Depth_range"].size() != 2))
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"Depth_range\"���ݻ���\"Depth_range\"���ݴ���\n", i + 1);
				exit(0);
			}
			double Depth_range[2];
			for (int j = 0; j < value["cameras"][i]["Depth_range"].size(); j++)
			{
				Depth_range[j] = value["cameras"][i]["Depth_range"][j].asDouble();
			}
			Rnears.push_back(Depth_range[0]);
			Rfars.push_back(Depth_range[1]);

			if (value["cameras"][i]["BitDepthColor"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"BitDepthColor\"���ݣ�\n", i + 1);
				exit(0);
			}
			texturebitDepths.push_back(value["cameras"][i]["BitDepthColor"].asUInt());

			if (value["cameras"][i]["BitDepthDepth"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"BitDepthDepth\"���ݣ�\n", i + 1);
				exit(0);
			}
			depthbitDepths.push_back(value["cameras"][i]["BitDepthDepth"].asUInt());

			if (value["cameras"][i]["ColorSpace"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"ColorSpace\"���ݣ�\n", i + 1);
				exit(0);
			}
			if (value["cameras"][i]["ColorSpace"].asString() == "YUV420" || value["cameras"][i]["ColorSpace"].asString() == "yuv420")
			{
				textureYUVFormats.push_back(YUV420);
			}
			else
			{
				printf("����:��֧�ֵ�����ͼ��ɫ�ռ䣺%s��", value["cameras"][i]["ColorSpace"].asString().c_str());
				exit(0);
			}

			if (value["cameras"][i]["DepthColorSpace"].isNull())
			{
				printf("����:json�ļ��е�%d���ӵ�ȱ��\"DepthColorSpace\"���ݣ�\n", i + 1);
				exit(0);
			}
			if (value["cameras"][i]["DepthColorSpace"].asString() == "YUV420" || value["cameras"][i]["DepthColorSpace"].asString() == "yuv420")
			{
				depthYUVFormats.push_back(YUV420);
			}
			else
			{
				printf("����:��֧�ֵ����ͼ��ɫ�ռ䣺%s��", value["cameras"][i]["DepthColorSpace"].asString().c_str());
				exit(0);
			}

			if (projections[i] == "Equirectangular")
			{
				if (value["cameras"][i]["Hor_range"].isNull() || (value["cameras"][i]["Hor_range"].size() != 2))
				{
					printf("����:json�ļ��е�%d���ӵ�ȱ��\"Hor_range\"���ݻ���\"Hor_range\"���ݴ���\n", i + 1);
					exit(0);
				}
				double Hor_range[2];
				for (int j = 0; j < value["cameras"][i]["Hor_range"].size(); j++)
				{
					Hor_range[j] = value["cameras"][i]["Hor_range"][j].asDouble();
				}
				horRanges.push_back(std::make_pair(Hor_range[0], Hor_range[1]));

				if (value["cameras"][i]["Ver_range"].isNull() || (value["cameras"][i]["Ver_range"].size() != 2))
				{
					printf("����:json�ļ��е�%d���ӵ�ȱ��\"Ver_range\"���ݻ���\"Ver_range\"���ݴ���\n", i + 1);
					exit(0);
				}
				double Ver_range[2];
				for (int j = 0; j < value["cameras"][i]["Ver_range"].size(); j++)
				{
					Ver_range[j] = value["cameras"][i]["Ver_range"][j].asDouble();
				}
				verRanges.push_back(std::make_pair(Ver_range[0], Ver_range[1]));
			}
			else if (projections[i] == "Perspective")
			{
				if (value["cameras"][i]["Principle_point"].isNull() || (value["cameras"][i]["Principle_point"].size() != 2))
				{
					printf("����:json�ļ��е�%d���ӵ�ȱ��\"Principle_point\"���ݻ���\"Principle_point\"���ݴ���\n", i + 1);
					exit(0);
				}
				double Principle_point[2];
				for (int j = 0; j < value["cameras"][i]["Principle_point"].size(); j++)
				{
					Principle_point[j] = value["cameras"][i]["Principle_point"][j].asDouble();
				}
				deltaXs.push_back(Principle_point[0]);
				deltaYs.push_back(Principle_point[1]);


				if (value["cameras"][i]["Focal"].isNull() || (value["cameras"][i]["Focal"].size() != 2))
				{
					printf("����:json�ļ��е�%d���ӵ�ȱ��\"Focal\"���ݻ���\"Focal\"���ݴ���\n", i + 1);
					exit(0);
				}
				double Focal[2];
				for (int j = 0; j < value["cameras"][i]["Focal"].size(); j++)
				{
					Focal[j] = value["cameras"][i]["Focal"][j].asDouble();
				}
				fxs.push_back(Focal[0]);
				fys.push_back(Focal[1]);
			}
			else
			{
				printf("����:��֧�ֵ�ͶӰ��ʽ��%s��\n", projections[i].c_str());
				printf("ֻ֧��һ��ͶӰ��ʽ��\n");
				printf("\tEquirectangular\n\tPerspective\n");
				exit(0);
			}
		}
		ifp.close();
	}
	void PrintParas()
	{
		printf("\n������\n");
		printf("�����ļ�����·����%s\n", (outputPlyFolder+"\\").c_str());
		printf("json�ļ�·����%s\n", jsonFilePath.c_str());
	}
	void PrintHelp()
	{
		std::cout<< "\n--����" << std::endl;
		printf("DT2PC.exe [-o Ŀ¼]  -j ·��\n\n");
		printf("\t-o\t\t�����ļ���Ҫ�����Ŀ¼��\n");
		printf("\t-j\t\tjson�ļ�����·��\n");
	}
};


class YUV
{
private:
	std::vector< std::vector <uint16_t> > *Y;
	std::vector< std::vector <uint16_t> > *U;
	std::vector< std::vector <uint16_t> > *V;
	uint16_t W;
	uint16_t H;
	uint8_t bitDepth;
	YUVformat cs;
public:
	std::vector< std::vector <uint16_t> > &operator[](size_t i)
	{
		assert(i <= 2);
		switch (i)
		{
		case 0: return *Y;
		case 1: return *U;
		case 2: return *V;
		}
	}
	YUV()
	{
		W = 0; H = 0; bitDepth = 0;
		Y = NULL; U = NULL; V = NULL;
	}
	YUV(uint16_t pW, uint16_t pH, uint8_t pbitDepth, YUVformat pcs):W(pW),H(pH),bitDepth(pbitDepth),cs(pcs)
	{
		assert(bitDepth <= 16);
		uint16_t emptyVal= 1 << (bitDepth - 1);
		Y = new std::vector< std::vector <uint16_t> >(H, std::vector <uint16_t>(W, emptyVal));
		U = new std::vector< std::vector <uint16_t> >(H, std::vector <uint16_t>(W, emptyVal));
		V = new std::vector< std::vector <uint16_t> >(H, std::vector <uint16_t>(W, emptyVal));
	}
	~YUV()
	{
		if (Y) delete Y;
		if (U) delete U;
		if (V) delete V;
	}
	
	void YUVRead(std::string fileName,size_t frameCounter)
	{
		assert(Y&&U&&V);
		if (cs == YUV420)
		{
			uint16_t *FileBuffer = (uint16_t*)malloc(H*W * sizeof(uint16_t) * 3 / 2);
			long long offset = H*W * 3 * frameCounter;
			
			FILE *Fp = fopen(fileName.c_str(), "rb");
			_fseeki64(Fp, offset, SEEK_SET);
			if (fread(FileBuffer, W*H * 3, 1, Fp) != 1)
			{
				printf("can't open the YUV File %s \n", fileName.c_str());
				exit(1);
			}
			fclose(Fp);

			for (int j = 0; j < H; j++)
			{
				for (int i = 0; i < W; i++)
				{
					int address_y = j * W + i;
					int address_u = H*W + (int)(j / 2)*(W / 2) + (int)(i / 2);
					int address_v = H*W * 5 / 4 + (int)(j / 2)*(W / 2) + (int)(i / 2);

					(*Y)[j][i] = FileBuffer[address_y];
					(*U)[j][i] = FileBuffer[address_u];
					(*V)[j][i] = FileBuffer[address_v];
					
				}
			}
			free(FileBuffer);
			return;
		}
		else
		{
			std::cout<<"Haven't relized reading such format YUV File!" << std::endl;
			exit(1);
		}
	}
	void YUVWrite(std::string fileName)
	{
		assert(Y&&U&&V);
		FILE *Fp = fopen(fileName.c_str(), "ab");

		for (int i = 0; i <H; i++)
		{
			for (int j = 0; j < W; j++)
			{
				fwrite(&((*Y)[i][j]), sizeof(uint16_t), 1, Fp);
			}
		}
		for (int i = 0; i < H; i += 2)
		{
			for (int j = 0; j < W; j += 2)
			{
				fwrite(&((*U)[i][j]), sizeof(uint16_t), 1, Fp);
			}
		}
		for (int i = 0; i < H; i += 2)
		{
			for (int j = 0; j < W; j += 2)
			{
				fwrite(&((*V)[i][j]), sizeof(uint16_t), 1, Fp);
			}
		}
		fclose(Fp);

		return;

	}
	bool YUV444To420(){return true;}
	bool YUV420To444(){return true;}
	void SetY(uint16_t m, uint16_t n, uint16_t val){assert(m < H&&n < W);(*Y)[m][n] = val;}
	void SetU(uint16_t m, uint16_t n, uint16_t val){assert(m < H&&n < W);(*U)[m][n] = val;}
	void SetV(uint16_t m, uint16_t n, uint16_t val){assert(m < H&&n < W);(*V)[m][n] = val;}
	uint16_t GetY(uint16_t m, uint16_t n) {return (*Y)[m][n];}
	uint16_t GetU(uint16_t m, uint16_t n) {return (*U)[m][n];}
	uint16_t GetV(uint16_t m, uint16_t n) {return (*V)[m][n];}
	
};

class CErpFrame
{
private:
	YUV *texture;
	YUV *depth;
	//������Ϣ
	std::string viewName;
	uint16_t width;
	uint16_t height;
	//����ͼ��Ϣ
	uint16_t texturebitDepth;
	YUVformat textureYUVFormat;
	//���ͼ��Ϣ
	double Rnear;
	double Rfar;
	YUVformat depthYUVFormat;
	uint16_t depthbitDepth;
	long vmax;

	Vector3<double> T;
	Vector3<double> R;

	double deltaX;
	double deltaY;
	double fx;
	double fy;


public:
	CErpFrame()
	{
		texture = NULL;
		depth = NULL;
		T = Vector3<double>(0.0, 0.0, 0.0);
		R = Vector3<double>(0.0, 0.0, 0.0);
	}
	CErpFrame(Vector3<double> pT, Vector3<double> pR) :T(pT), R(pR) { texture = NULL; depth = NULL; }
	CErpFrame(std::string pviewName, uint16_t pwidth, uint16_t pheight,YUVformat pdepthYUVFormat, uint8_t pdepthbitDepth,
		YUVformat ptextureYUVFormat, uint8_t ptexturebitDepth)
	{
		viewName = pviewName;
		width = pwidth;
		height = pheight;
		//����ͼ��Ϣ
		texturebitDepth = ptexturebitDepth;
		textureYUVFormat = ptextureYUVFormat;
		//���ͼ��Ϣ
		depthYUVFormat = pdepthYUVFormat;
		depthbitDepth = pdepthbitDepth;

		texture = new YUV(width, height, texturebitDepth, textureYUVFormat);
		depth = new YUV(width, height, depthbitDepth, depthYUVFormat);
	}
	CErpFrame(std::string pviewName, uint16_t pwidth, uint16_t pheight, double pRnear, double pRfar, YUVformat pdepthYUVFormat, uint8_t pdepthbitDepth,
		YUVformat ptextureYUVFormat, uint8_t ptexturebitDepth, Vector3<double> pT, Vector3<double> pR,double pdeltaX,double pdeltaY,
		double pfx,double pfy)
	{
		viewName = pviewName;
		width = pwidth;
		height = pheight;
		//����ͼ��Ϣ
		texturebitDepth = ptexturebitDepth;
		textureYUVFormat = ptextureYUVFormat;
		//���ͼ��Ϣ
		Rnear = pRnear;
		Rfar = pRfar;
		depthYUVFormat = pdepthYUVFormat;
		depthbitDepth = pdepthbitDepth;
		T = pT;
		R = pR;
		deltaX = pdeltaX;
		deltaY = pdeltaY;
		fx = pfx;
		fy = pfy;
		if (texture) delete texture;
		if (depth) delete depth;
		texture = new YUV(width, height, texturebitDepth, textureYUVFormat);
		depth = new YUV(width, height, depthbitDepth, depthYUVFormat);
		vmax = 1 << depthbitDepth - 1;
	}
	~CErpFrame()
	{
		if (texture) delete texture;
		if (depth) delete depth;
	}

	YUV *getTexture() { return texture; }
	YUV *getDepth() { return depth; }
	std::string &getViewName() { return viewName; }
	uint16_t &getWidth() { return width; }
	uint16_t &getHeight() { return height; }
	uint16_t &getTexturebitDepth() { return texturebitDepth; }
	YUVformat &getTextureYUVFormat() { return textureYUVFormat; }
	double &getRnear() { return Rnear; }
	double &getRfar() { return Rfar; }
	YUVformat &getDepthYUVFormat() { return depthYUVFormat; }
	uint16_t &getDepthbitDepth() { return depthbitDepth; }
	long &getVmax() { return vmax; }
	Vector3<double> &getT() { return T; }
	Vector3<double> &getR() { return R; }
	double &getdeltaX() { return deltaX; }
	double &getdeltaY() { return deltaY; }
	double &getfx() { return fx; }
	double &getfy() { return fy; }

	void setTAndR(Vector3<double> pT, Vector3<double> pR) { T = pT; R = pR; }
	void textureYUV420ToYUV444()
	{
		assert(texture&&depth&&textureYUVFormat == YUV420);
		texture->YUV420To444();
	}
	void textureYUV444ToYUV420()
	{
		assert(texture&&depth&&textureYUVFormat == YUV444);
		texture->YUV444To420();
	}
	void read(std::string textureFileName, std::string depthFileName, size_t frameCounter)
	{
		if (texture) delete texture;
		if (depth) delete depth;
		texture = new YUV(width, height, texturebitDepth, textureYUVFormat);
		depth = new YUV(width, height, depthbitDepth, depthYUVFormat);
		texture->YUVRead(textureFileName, frameCounter);
		depth->YUVRead(depthFileName, frameCounter);
		return;
	}
};

int sign(double Z)
{
	return (Z > 0 ? 1 : -1);
}

#endif // !erp_h

