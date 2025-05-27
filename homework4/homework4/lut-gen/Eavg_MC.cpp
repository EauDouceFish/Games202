#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>
#include <random>
#include "vec.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int resolution = 128;
int channel = 3;

typedef struct samplePoints {
    std::vector<Vec3f> directions;
	std::vector<float> PDFs;
}samplePoints;

void setRGB(int x, int y, float alpha, unsigned char *data) {
	data[3 * (resolution * x + y) + 0] = uint8_t(alpha);
    data[3 * (resolution * x + y) + 1] = uint8_t(alpha);
    data[3 * (resolution * x + y) + 2] = uint8_t(alpha);
}

void setRGB(int x, int y, Vec3f alpha, unsigned char *data) {
	data[3 * (resolution * x + y) + 0] = uint8_t(alpha.x);
    data[3 * (resolution * x + y) + 1] = uint8_t(alpha.y);
    data[3 * (resolution * x + y) + 2] = uint8_t(alpha.z);
}

samplePoints squareToCosineHemisphere(int sample_count){
    samplePoints samlpeList;
    const int sample_side = static_cast<int>(floor(sqrt(sample_count)));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rng(0.0, 1.0);
    for (int t = 0; t < sample_side; t++) {
        for (int p = 0; p < sample_side; p++) {
            double samplex = (t + rng(gen)) / sample_side;
            double sampley = (p + rng(gen)) / sample_side;
            
            double theta = 0.5f * acos(1 - 2*samplex);
            double phi =  2 * PI * sampley;
            Vec3f wi = Vec3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
            float pdf = wi.z / PI;
            
            samlpeList.directions.push_back(wi);
            samlpeList.PDFs.push_back(pdf);
        }
    }
    return samlpeList;
}

Vec3f getEmu(int x, int y, int alpha, unsigned char *data, float NdotV, float roughness) {
    return Vec3f(data[3 * (resolution * x + y) + 0],
                 data[3 * (resolution * x + y) + 1],
                 data[3 * (resolution * x + y) + 2]);
}

Vec3f IntegrateEmu(Vec3f V, float roughness, float NdotV, Vec3f Ei) {
    Vec3f Eavg = Vec3f(0.0f);
    const int sample_count = 1024;
    // Vec3f N = Vec3f(0.0, 0.0, 1.0);

    // samplePoints sampleList = squareToCosineHemisphere(sample_count);
    for (int i = 0; i < sample_count; i++) {
        // 本质上是mu
        // Vec3f L = sampleList.directions[i];
        // Vec3f H = normalize(V + L);

        // float NoL = std::max(L.z, 0.0f);
        // float NoH = std::max(H.z, 0.0f);
        // float VoH = std::max(dot(V, H), 0.0f);
        // float NoV = std::max(dot(N, V), 0.0f);

        // TODO: To calculate Eavg here
        Vec3f E_mu = Ei;
        float mu = NdotV;
        Eavg += E_mu * mu * 2.0f;
    }
    // 除以行宽
    return Eavg / sample_count;
}

void printProgress(int current, int total) {
    int barWidth = 50; // 进度条宽度
    float progress = float(current) / float(total);
    int pos = int(barWidth * progress);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %";
    std::cout.flush();
}

int main() {
    unsigned char *Edata = stbi_load("./GGX_E_MC_LUT.png", &resolution, &resolution, &channel, 3);
    if (Edata == NULL) 
    {
		std::cout << "ERROE_FILE_NOT_LOAD" << std::endl;
		return -1;
	}
	else 
    {
		std::cout << "Resolution: " << resolution << "*" << resolution << " Channel: " << channel << std::endl;
        std::cout << "Creating the LUT..." << std::endl;
        // | -----> mu(j)
        // | 
        // | rough（i）
        // flip it if you want to write the data on picture 
        uint8_t* data = new uint8_t[resolution * resolution * 3];
        float step = 1.0 / resolution;
        Vec3f Eavg = Vec3f(0.0);
		for (int i = 0; i < resolution; i++) 
        {
            float roughness = step * (static_cast<float>(i) + 0.5f);
			for (int j = 0; j < resolution; j++) 
            {
                float NdotV = step * (static_cast<float>(j) + 0.5f);
                Vec3f V = Vec3f(std::sqrt(1.f - NdotV * NdotV), 0.f, NdotV);

                Vec3f Ei = getEmu((resolution - 1 - i), j, 0, Edata, NdotV, roughness);
                Eavg += IntegrateEmu(V, roughness, NdotV, Ei) * step;
                setRGB(i, j, 0.0, data);
			}

            for(int k = 0; k < resolution; k++)
            {
                setRGB(i, k, Eavg, data);
            }

            Eavg = Vec3f(0.0);
            printProgress(i + 1, resolution);
		}

		stbi_flip_vertically_on_write(true);
		stbi_write_png("GGX_Eavg_MC_LUT.png", resolution, resolution, channel, data, 0);
	}
	stbi_image_free(Edata);
    return 0;
}