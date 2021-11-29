#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <geometry.hh>

#define UseBackground 0

namespace PC
{
    const unsigned W = 1920, H = 1080;
    std::vector<Vec3f> FrameBuffer(W*H);
    std::string OutFilename;
    
    const unsigned FOV = M_PI/2.;
    const unsigned ReflectionDepth = 3;
    
    int BackgroundWidth, BackgroundHeight;
    std::vector<Vec3f> Background;
}

class LightSource
{
    public:
        LightSource(const Vec3f &Position, const Vec3f &Color, const float &Intensity) :
        Position(Position), Color(Color), Intensity(Intensity) {}
    
        Vec3f Position;
        Vec3f Color;
        float Intensity;
};

class Material
{
    public:
        Material(const float &RefractiveIndex, const Vec4f &Albedo, const Vec3f &Color, const float &SpecularExpoent) :
        RefractiveIndex(RefractiveIndex), Albedo(Albedo), DiffuseColor(Color), SpecularExpoent(SpecularExpoent){}
        
        Material() : RefractiveIndex(1), Albedo(1,0,0,0), DiffuseColor(), SpecularExpoent() {}
    
        float RefractiveIndex;
        Vec4f Albedo;
        Vec3f DiffuseColor;
        float SpecularExpoent;
};

class Sphere
{
    public:
        Sphere(const Vec3f &Center, const float &Radius,
               const Material &SphereMaterial) : Center(Center), Radius(Radius),
                                                 SphereMaterial(SphereMaterial) {}
        
        bool
        RayIntersect(const Vec3f &Src, const Vec3f &Dir, float &t0) const {
            Vec3f L = Center - Src;
            float TCA = L*Dir;
            float D2 = L*L-TCA*TCA;
            if (D2 > Radius*Radius) return false;
            
            float THC = sqrtf(Radius*Radius-D2);
            t0 		 = TCA-THC;
            float t1 = TCA+THC;
            if (t0<0) t0 = t1;
            if (t0<0) return false;
            
            return true;
        }
        
        Vec3f Center;
        float Radius;
        Material SphereMaterial;
};

Vec3f
Refract(const Vec3f &I, const Vec3f &N, const float &RefractiveIndex)
{
    float CosI = - std::max(-1.f, std::min(1.f, I*N));
    float EtaI = 1, EtaT = RefractiveIndex;
    Vec3f n = N;
    if (CosI<0) {
        CosI = -CosI;
        std::swap(EtaI, EtaT);
        n = -N;
    }
    float Eta = EtaI/EtaT;
    float K = 1-Eta*Eta*(1-CosI*CosI);
    return K<0?Vec3f(0,0,0):I*Eta+n*(Eta*CosI-sqrtf(K));
}

Vec3f
Reflect(const Vec3f &I, const Vec3f &N)
{
    return I-N*2.f*(I*N);
}
 
bool
SceneIntersect(const Vec3f &Src, const Vec3f &Dir, const std::vector<Sphere> &Spheres,
               Vec3f &Hit, Vec3f &N, Material &SphereMaterial)
{
    float SpheresDistance = std::numeric_limits<float>::max();
    
    for (auto &Sphere : Spheres) {
        float Distance;
        if (Sphere.RayIntersect(Src, Dir, Distance) && Distance < SpheresDistance) {
            SpheresDistance = Distance;
            Hit = Src+Dir*Distance;
            N = (Hit-Sphere.Center).normalize();
            SphereMaterial = Sphere.SphereMaterial;
        }
    }
    
    return SpheresDistance<3000;
}

Vec3f
CastRay(const Vec3f &Src, const Vec3f &Dir, const std::vector<Sphere> &Spheres,
        const std::vector<LightSource> &LightSources, unsigned Depth=0)
{
    Vec3f Point, N;
    Material SphereMaterial;
    
    if (Depth>PC::ReflectionDepth || !SceneIntersect(Src, Dir, Spheres, Point, N, SphereMaterial)) {
        if (UseBackground) {
            int a = std::max(0, std::min(PC::BackgroundWidth -1, static_cast<int>((atan2(Dir.z, Dir.x)/(2*M_PI)+.5)*PC::BackgroundWidth)));
            int b = std::max(0, std::min(PC::BackgroundHeight-1, static_cast<int>(acos(Dir.y)/M_PI*PC::BackgroundHeight)));
            return PC::Background[a+b*PC::BackgroundWidth];
        } else {
            return Vec3f(0.0, 0.0, 0.0);
        }
    }
    
    Vec3f ReflectDirection = Reflect(Dir, N).normalize();
    Vec3f RefractDirection = Refract(Dir, N, SphereMaterial.RefractiveIndex).normalize();
    Vec3f RefractSource = RefractDirection*N<0?Point-N*1e-3:Point+N*1e-3;
    Vec3f ReflectSource = ReflectDirection*N<0?Point-N*1e-3:Point+N*1e-3;
    Vec3f ReflectColor = CastRay(ReflectSource,  ReflectDirection, Spheres, LightSources, Depth+1);
    Vec3f RefractColor = CastRay(RefractSource, RefractDirection, Spheres, LightSources, Depth+1);
    
    float DiffuseLightIntensity = 0, SpecularLightIntensity = 0;
    for (auto &LightSource : LightSources) {
        Vec3f LightDirection = (LightSource.Position - Point).normalize();
        float LightDistance  = (LightSource.Position - Point).norm();
        
        Vec3f ShadowSource = LightDirection*N<0?Point-N*1e-3:Point+N*1e-3; // Checks if the point is in shadow of LightSource
        Vec3f ShadowPoint, ShadowN;
        Material TemporaryMaterial;
        if ((SceneIntersect(ShadowSource, LightDirection, Spheres, ShadowPoint, ShadowN, TemporaryMaterial)) &&
            (ShadowPoint-ShadowSource).norm()<LightDistance) continue;
        
        DiffuseLightIntensity += LightSource.Intensity*std::max(0.f, LightDirection*N);
        SpecularLightIntensity += powf(std::max(0.f, Reflect(LightDirection, N)*Dir),
                                  SphereMaterial.SpecularExpoent) * LightSource.Intensity;
    }
    
    return SphereMaterial.DiffuseColor * DiffuseLightIntensity * SphereMaterial.Albedo[0] +
           Vec3f(1., 1., 1.) * SpecularLightIntensity * SphereMaterial.Albedo[1] +
           ReflectColor*SphereMaterial.Albedo[2] + RefractColor*SphereMaterial.Albedo[3];
}

void
Render(const std::vector<Sphere> &Spheres, const std::vector<LightSource> &LightSources)
{
    #pragma omp parallel for
    for (unsigned ScreenY = 0; ScreenY < PC::H; ScreenY++)
        for (unsigned ScreenX = 0; ScreenX < PC::W; ScreenX++) {
            float X =  (2*(ScreenX+0.5)/(float)PC::W-1)*tan(PC::FOV/2.)*PC::W/(float)PC::H;
            float Y = -(2*(ScreenY+0.5)/(float)PC::H-1)*tan(PC::FOV/2.);
            Vec3f Dir = Vec3f(X, Y, -1).normalize();
            
            PC::FrameBuffer[ScreenX+ScreenY*PC::W] = CastRay(Vec3f(0,0,0), Dir, Spheres, LightSources);
        }
    
    // Writing RGB pixels to the image file
    std::vector<unsigned char> PixMap(PC::W*PC::H*3);
    for (unsigned x = 0; x < PC::H*PC::W; ++x) {
        Vec3f &C = PC::FrameBuffer[x];
        float Max = std::max(C[0], std::max(C[1], C[2]));
        if (Max>1) C = C*(1./Max);
        for (unsigned y = 0; y<3; y++)
            PixMap[x*3+y] = (unsigned char)(255*std::max(0.f, std::min(1.f, PC::FrameBuffer[x][y])));
    }
    stbi_write_jpg(PC::OutFilename.c_str(), PC::W, PC::H, 3, PixMap.data(), 100);
}

int
main()
{
    // Background image rendering stuff
    if (UseBackground) {
        int n = -1;
        unsigned char *PixMap = stbi_load("./res/background.jpg", &PC::BackgroundWidth, &PC::BackgroundHeight, &n, 0);
        if (!PixMap || 3!=n) {
            std::cerr << "Error: couldn't load the background map\n";
            exit(-1);
        }
        
        PC::Background = std::vector<Vec3f>(PC::BackgroundWidth*PC::BackgroundHeight);
        for (int y=PC::BackgroundHeight-1; y>=0; y--)
            for (int x=0; x<PC::BackgroundWidth; x++)
                PC::Background[x+y*PC::BackgroundWidth] = Vec3f(PixMap[(x+y*PC::BackgroundWidth)*3+0],
                                                                PixMap[(x+y*PC::BackgroundWidth)*3+1],
                                                                PixMap[(x+y*PC::BackgroundWidth)*3+2])*(1/255.);
        stbi_image_free(PixMap);
    }

    // The colors are in a range of 0 to 1, instead of 0 to 255.
    // On GIMP, you can convert the 255-range RGB values to 100-range
    // RGB values and divide them by 100.
    Material  WhitePlastic(0.5, Vec4f(0.9, 0.5,  0.3, 0.0), Vec3f(1.0, 1.0, 1.0), 50.f);
    Material 	RedPlastic(0.5, Vec4f(0.9, 0.5,  0.3, 0.0), Vec3f(0.8, 0.1, 0.2), 50.f);
    Material      BlueWall(1.0, Vec4f(0.9, 0.1,  0.1, 0.0), Vec3f(0.2, 0.2, 1.0), 50.f);
    Material YellowPlastic(0.5, Vec4f(0.9, 0.5,  0.3, 0.0), Vec3f(1.0, 0.5, 0.2), 50.f);
    Material   PinkPlastic(0.5, Vec4f(0.9, 0.5,  0.3, 0.0), Vec3f(0.9, 0.1, 0.9), 50.f);
    Material 	    Mirror(1.0, Vec4f(0.0, 10.0, 0.6, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.f);
    
    std::vector<LightSource> LightSources = {
        LightSource(Vec3f( 25,  0,   185), Vec3f(1.0, 1.0, 1.0), 0.2),
        LightSource(Vec3f( 0,   0,   -50), Vec3f(1.0, 1.0, 1.0), 0.2),
        LightSource(Vec3f( 0,   0,    10), Vec3f(1.0, 1.0, 1.0), 0.3),
    };
    
    float Z[14] = {
        -20,  -30,  -60,
        -70,  -80,  -110,
        -120, -130, -160,
        -170, -190, -220,
    };
    
    float WallsZ = -50;
    float FrontWallZ = -2998;
    
    float Move = 0.1;
    
    for (unsigned Frame = 1; Frame < 1000; Frame++) {
        PC::OutFilename = "./frames/frame"+std::to_string(Frame)+".jpg";
        
        std::vector<Sphere> Spheres = {
            Sphere(Vec3f(-2, -3,  Z[1]), 3, RedPlastic),
            Sphere(Vec3f( 23, 10, Z[2]), 4, PinkPlastic),
            Sphere(Vec3f( 0,  5,  Z[1]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[0]), 2, RedPlastic),
            Sphere(Vec3f( 7,  5,  Z[1]), 3, Mirror),
            Sphere(Vec3f(-5, -1,  Z[2]), 2, Mirror),
            Sphere(Vec3f( 0,  5,  Z[3]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[3]), 2, YellowPlastic),
            Sphere(Vec3f(-2, -3,  Z[5]), 3, RedPlastic),
            Sphere(Vec3f( 23, 10, Z[6]), 4, PinkPlastic),
            Sphere(Vec3f( 0,  5,  Z[5]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[4]), 2, RedPlastic),
            Sphere(Vec3f( 7,  5,  Z[5]), 3, Mirror),
            Sphere(Vec3f(-5, -1,  Z[5]), 2, Mirror),
            Sphere(Vec3f( 0,  5,  Z[6]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[6]), 2, YellowPlastic),
            Sphere(Vec3f(-5, -1,  Z[8]), 2, Mirror),
            Sphere(Vec3f( 0,  5,  Z[9]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[9]), 2, YellowPlastic),
            Sphere(Vec3f(-2, -3,  Z[11]),3, RedPlastic),
            Sphere(Vec3f( 23, 10, Z[10]),4, PinkPlastic),
            Sphere(Vec3f( 0,  5,  Z[9]), 2, YellowPlastic),
            Sphere(Vec3f( 4, -2,  Z[8]), 2, RedPlastic),
            
            // Walls
            Sphere(Vec3f(-2998, 0, WallsZ),  2890, BlueWall), // Left
            Sphere(Vec3f( 2998, 0, WallsZ),  2890, BlueWall), // Right
            Sphere(Vec3f( 0, 2998, WallsZ),  2890, BlueWall), // Top
            Sphere(Vec3f( 0,-2998, WallsZ),  2890, BlueWall), // Bottom
            Sphere(Vec3f( 0, 0, FrontWallZ), 2250, BlueWall), // Front
        };
        
        std::cout << "Rendering frame " << Frame << "...\n";
        Render(Spheres, LightSources);
        
        for (auto &ZElement : Z)
            ZElement += Move;
        
        WallsZ +=     Move*5;
        FrontWallZ += Move*5;
    }
    
    return 0;
}
