#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>
#include "vec.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int resolution = 128;
int channel = 3;

Vec2f Hammersley(uint32_t i, uint32_t N) {
    uint32_t bits = (i << 16u) | (i >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    float rdi = float(bits) * 2.3283064365386963e-10;
    return {float(i) / float(N), rdi};
}

Vec3f ImportanceSampleGGX(Vec2f xi, Vec3f N, float roughness) {

    float a = roughness * roughness;

    // TODO: Copy the code from your previous work - Bonus 1

    float theta_m = atan(a * sqrt(xi.x) / sqrt(1 - xi.x));
    float phi_h = 2.0 * PI * xi.y;
    
    //TODO: from spherical space to cartesian space - Bonus 1
    Vec3f H_local;
    H_local.x = sin(theta_m) * cos(phi_h);
    H_local.y = sin(theta_m) * sin(phi_h);
    H_local.z = cos(theta_m);

    //TODO: tangent coordinates - Bonus 1
    Vec3f up        = abs(N.z) < 0.999 ? Vec3f(0.0, 0.0, 1.0) : Vec3f(1.0, 0.0, 0.0);
    Vec3f tangent   = normalize(cross(up, N));
    Vec3f bitangent = cross(N, tangent);

    //TODO: transform H from tangent space to world - Bonus 1
    Vec3f sample_vec = tangent * H_local.x + bitangent * H_local.y + N * H_local.z;
    
    return Vec3f(sample_vec);
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

Vec3f getEmu(int x, int y, int alpha, unsigned char *data, float NdotV, float roughness) {
    return Vec3f(data[3 * (resolution * x + y) + 0],
                 data[3 * (resolution * x + y) + 1],
                 data[3 * (resolution * x + y) + 2]);
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
    unsigned char *Edata = stbi_load("./GGX_E_LUT.png", &resolution, &resolution, &channel, 3);
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
        // Flip it, if you want the data written to the texture
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
            printProgress(i + 1, resolution);
            Eavg = Vec3f(0.0);
		}
		stbi_flip_vertically_on_write(true);
		stbi_write_png("GGX_Eavg_LUT.png", resolution, resolution, channel, data, 0);
	}
	stbi_image_free(Edata);
    return 0;
}