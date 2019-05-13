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

	//view 参数
	std::vector<std::string> viewNames;
	std::vector<uint16_t> widths;
	std::vector<uint16_t> heights;
	std::vector<Vector3<double>> Ts; //各个视点平移向量，由json文件解析得到（解析参数）
	std::vector<Vector3<double>> Rs; //各个视点旋转矩阵，由json文件文件解析得到（解析参数）
	std::vector<double> Rnears;
	std::vector<double> Rfars;
	std::vector<std::string> projections;

	//纹理图信息
	std::vector<std::string> textureFileNames;
	std::vector<uint16_t> texturebitDepths;
	std::vector<YUVformat> textureYUVFormats;
	//深度图信息
	std::vector<std::string> depthFileNames;
	std::vector<YUVformat> depthYUVFormats;
	std::vector<uint16_t> depthbitDepths;

	//全景图信息
	std::vector<std::pair<double, double>> horRanges;
	std::vector<std::pair<double, double>> verRanges;
	//透视图信息
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
			else { printf("警告：错误的的输入参数：%s!\n", argv[iArgIdx]);  PrintHelp();}
		}
		//检查命令行参数
		if (outputPlyFolder.empty()) { printf("警告：未输入点云文件保存路径!\n"); outputPlyFolder = "."; }
		if (jsonFilePath.empty()) { printf("错误：未输入json文件路径!\n"); PrintHelp(); exit(0); }

		//解析json文件
		Json::Value value;
		Json::Reader jsonReader;
		std::ifstream ifp;
		ifp.open(jsonFilePath);
		if (!ifp.is_open())
		{
			printf("错误：json文件打开失败！\n");
			exit(0);
		}

		jsonReader.parse(ifp, value);
		if (value["Content_name"].isNull())
		{
			printf("错误:json文件中缺少\"Content_name\"数据！\n");
			exit(0);
		}
		ContentName = value["Content_name"].asString();

		if (value["Frames_number"].isNull())
		{
			printf("错误:json文件中缺少\"Frames_number\"数据！\n");
			exit(0);
		}
		frameNumber = value["Frames_number"].asUInt64();

		if (value["cameras"].isNull())
		{
			printf("错误:json文件中缺少\"cameras\"数据！\n");
			exit(0);
		}
		for (int i = 0; i < value["cameras"].size(); i++)
		{
			if (value["cameras"][i]["Name"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"Name\"数据！\n", i + 1);
				exit(0);
			}
			viewNames.push_back(value["cameras"][i]["Name"].asString());

			if (value["cameras"][i]["textureFilePath"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"textureFilePath\"数据！\n", i + 1);
				exit(0);
			}
			textureFileNames.push_back(value["cameras"][i]["textureFilePath"].asString());

			if (value["cameras"][i]["depthFilePath"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"depthFilePath\"数据！\n", i + 1);
				exit(0);
			}
			depthFileNames.push_back(value["cameras"][i]["depthFilePath"].asString());

			if (value["cameras"][i]["Position"].isNull() || (value["cameras"][i]["Position"].size() != 3))
			{
				printf("错误:json文件中第%d个视点缺少\"Position\"数据或者\"Position\"数据错误！\n", i + 1);
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
				printf("错误:json文件中第%d个视点缺少\"Rotation\"数据或者\"Rotation\"数据错误！\n", i + 1);
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
				printf("错误:json文件中第%d个视点缺少\"Resolution\"数据或者\"Resolution\"数据错误！\n", i + 1);
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
				printf("错误:json文件中第%d个视点缺少\"Projection\"数据！\n", i + 1);
				exit(0);
			}
			projections.push_back(value["cameras"][i]["Projection"].asString());

			if (value["cameras"][i]["Depth_range"].isNull() || (value["cameras"][i]["Depth_range"].size() != 2))
			{
				printf("错误:json文件中第%d个视点缺少\"Depth_range\"数据或者\"Depth_range\"数据错误！\n", i + 1);
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
				printf("错误:json文件中第%d个视点缺少\"BitDepthColor\"数据！\n", i + 1);
				exit(0);
			}
			texturebitDepths.push_back(value["cameras"][i]["BitDepthColor"].asUInt());

			if (value["cameras"][i]["BitDepthDepth"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"BitDepthDepth\"数据！\n", i + 1);
				exit(0);
			}
			depthbitDepths.push_back(value["cameras"][i]["BitDepthDepth"].asUInt());

			if (value["cameras"][i]["ColorSpace"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"ColorSpace\"数据！\n", i + 1);
				exit(0);
			}
			if (value["cameras"][i]["ColorSpace"].asString() == "YUV420" || value["cameras"][i]["ColorSpace"].asString() == "yuv420")
			{
				textureYUVFormats.push_back(YUV420);
			}
			else
			{
				printf("错误:不支持的纹理图颜色空间：%s！", value["cameras"][i]["ColorSpace"].asString().c_str());
				exit(0);
			}

			if (value["cameras"][i]["DepthColorSpace"].isNull())
			{
				printf("错误:json文件中第%d个视点缺少\"DepthColorSpace\"数据！\n", i + 1);
				exit(0);
			}
			if (value["cameras"][i]["DepthColorSpace"].asString() == "YUV420" || value["cameras"][i]["DepthColorSpace"].asString() == "yuv420")
			{
				depthYUVFormats.push_back(YUV420);
			}
			else
			{
				printf("错误:不支持的深度图颜色空间：%s！", value["cameras"][i]["DepthColorSpace"].asString().c_str());
				exit(0);
			}

			if (projections[i] == "Equirectangular")
			{
				if (value["cameras"][i]["Hor_range"].isNull() || (value["cameras"][i]["Hor_range"].size() != 2))
				{
					printf("错误:json文件中第%d个视点缺少\"Hor_range\"数据或者\"Hor_range\"数据错误！\n", i + 1);
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
					printf("错误:json文件中第%d个视点缺少\"Ver_range\"数据或者\"Ver_range\"数据错误！\n", i + 1);
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
					printf("错误:json文件中第%d个视点缺少\"Principle_point\"数据或者\"Principle_point\"数据错误！\n", i + 1);
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
					printf("错误:json文件中第%d个视点缺少\"Focal\"数据或者\"Focal\"数据错误！\n", i + 1);
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
				printf("错误:不支持的投影格式：%s！\n", projections[i].c_str());
				printf("只支持一下投影格式：\n");
				printf("\tEquirectangular\n\tPerspective\n");
				exit(0);
			}
		}
		ifp.close();
	}
	void PrintParas()
	{
		printf("\n参数：\n");
		printf("点云文件保存路径：%s\n", (outputPlyFolder+"\\").c_str());
		printf("json文件路径：%s\n", jsonFilePath.c_str());
	}
	void PrintHelp()
	{
		std::cout<< "\n--帮助" << std::endl;
		printf("DT2PC.exe [-o 目录]  -j 路径\n\n");
		printf("\t-o\t\t点云文件将要保存的目录。\n");
		printf("\t-j\t\tjson文件所在路径\n");
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
	//共有信息
	std::string viewName;
	uint16_t width;
	uint16_t height;
	//纹理图信息
	uint16_t texturebitDepth;
	YUVformat textureYUVFormat;
	//深度图信息
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
		//纹理图信息
		texturebitDepth = ptexturebitDepth;
		textureYUVFormat = ptextureYUVFormat;
		//深度图信息
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
		//纹理图信息
		texturebitDepth = ptexturebitDepth;
		textureYUVFormat = ptextureYUVFormat;
		//深度图信息
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

