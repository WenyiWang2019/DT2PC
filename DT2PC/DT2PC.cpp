// DT2PC.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <string>
#include "common.h"
#include "PointSet.h"


void ToPointCloud(CErpFrame &view, PointSet3 &pointCloud, size_t viewsNum, Parameter &para);

int main(int argc, char *argv[])
{
	Parameter para;
	para.Parse(argc, argv);
	para.PrintParas();

	for (size_t i = 0; i < para.frameNumber; i++) //帧数
	{
		PointSet3 synPointCloud;
		printf("\n正在转换第%d帧......\n", i + 1);
		for (int j = 0; j < para.viewNames.size(); j++) //每帧视点数
		{
			printf("\t正在合并第%d个视点......\n", j + 1);

			CErpFrame view(para.viewNames[j], para.widths[j], para.heights[j], para.depthYUVFormats[j],
				para.depthbitDepths[j], para.textureYUVFormats[j], para.texturebitDepths[j]);
			view.read(para.textureFileNames[j], para.depthFileNames[j], i);
			view.textureYUV420ToYUV444();

			PointSet3 pointCloud;
			ToPointCloud(view, pointCloud, j, para);
			pointCloud.ToWorldCoordinates(para.Rs[j], para.Ts[j]);
			synPointCloud.addPoints(pointCloud);
		}

		synPointCloud.Voxelize(16);
		synPointCloud.sort();
		synPointCloud.RemoveRepeatePoints();

		if (para.texturebitDepths[0] != 8 && para.texturebitDepths[0] != 10)
		{
			printf("错误：暂时不支持该纹理精度：%d\n", para.texturebitDepths[0]);
			printf("\t提示：支持的纹理精度：8bit\t10bit%d\n");
			exit(0);
		}
		if(para.texturebitDepths[0]==10) synPointCloud.ConvertColorFrom10bTo8b();
		synPointCloud.convertYUVToRGB();

		char frameNum[8];
		sprintf(frameNum, "%04d", i);
		std::string synPointCloudName = para.outputPlyFolder + "\\syn_" + para.ContentName + "_F" + std::string(frameNum) + ".ply";
		printf("\t正在将点云文件写入磁盘......\n");
		synPointCloud.write(synPointCloudName, true);
	}

	return 0;
}

void ToPointCloud(CErpFrame &view, PointSet3 &pointCloud, size_t viewsNum, Parameter &para)
{
	auto W = para.widths[viewsNum];
	auto H = para.heights[viewsNum];
	long vmax = (1 << para.depthbitDepths[viewsNum]) - 1;
	double Rnear = para.Rnears[viewsNum];
	double Rfar = para.Rfars[viewsNum];
	Point3D T = para.Ts[viewsNum];
	Point3D Rotation = para.Rs[viewsNum];

	auto depth = view.getDepth();
	auto texture = view.getTexture();
	Point3D point3D;
	Color3B Color3B;
	uint16_t Depth;
	size_t piontCount = pointCloud.getPointCount();
	pointCloud.resize(piontCount + W * H);
		
	if (para.projections[viewsNum] == "Equirectangular")
	{
		double H_left = para.horRanges[viewsNum].first;
		double H_right = para.horRanges[viewsNum].second;
		double V_up = para.verRanges[viewsNum].first;
		double V_down = para.verRanges[viewsNum].second;

		double phi;
		double theta;
		double R;
		for (size_t m = 0; m < H; m++)
		{
			for (size_t n = 0; n < W; n++)
			{
				Depth = (*depth)[0][m][n];

				phi = ((n + 0.5) / W - 0.5) * (fabs(H_right - H_left) / 180) * PI;
				theta = (0.5 - (m + 0.5) / H)*(fabs(V_up - V_down) / 180)*PI;

				if (Depth != 0)
				{
					R = vmax * Rnear*Rfar / (Depth*(Rfar - Rnear) + vmax * Rnear);

					point3D[0] = R * cos(theta)*cos(phi);
					point3D[1] = -R * cos(theta)*sin(phi);
					point3D[2] = R * sin(theta);

					Color3B[0] = (*texture)[0][m][n];
					Color3B[1] = (*texture)[1][m][n];
					Color3B[2] = (*texture)[2][m][n];

					pointCloud.setPoint(piontCount, point3D);
					pointCloud.setColor(piontCount, Color3B);
					piontCount++;
				}
			}
		}
		pointCloud.resize(piontCount);
	}
	else if (para.projections[viewsNum] == "Perspective")
	{
		double deltaX = para.deltaXs[viewsNum];
		double deltaY = para.deltaYs[viewsNum];
		double fx = para.fxs[viewsNum];
		double fy = para.fys[viewsNum];

		for (size_t m = 0; m < H; m++)
		{
			for (size_t n = 0; n < W; n++)
			{
				Depth = (*depth)[0][m][n];

				if (Depth != 0)
				{
					//计算三维坐标
					point3D[0] = -vmax * Rnear*Rfar / (Depth*(Rfar - Rnear) + vmax * Rnear);
					point3D[1] = -point3D[0] * (n - deltaX) / fx;
					point3D[2] = point3D[0] * (m - deltaY) / fy;

					//获得纹理
					Color3B[0] = (*texture)[0][m][n];
					Color3B[1] = (*texture)[1][m][n];
					Color3B[2] = (*texture)[2][m][n];

					pointCloud.setPoint(piontCount, point3D);
					pointCloud.setColor(piontCount, Color3B);
					piontCount++;
				}
			}
		}
		pointCloud.resize(piontCount);
	}
	else
	{
		printf("错误:不支持的投影格式：%s！\n", para.projections[viewsNum].c_str());
		printf("只支持一下投影格式：\n");
		printf("\tEquirectangular\n\tPerspective\n");
		exit(0);
	}

	
}



