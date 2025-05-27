#ifdef GL_ES
precision mediump float;
#endif

// Phong related variables
uniform sampler2D uSampler;
uniform vec3 uKd;
uniform vec3 uKs;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;
uniform vec3 uLightIntensity;

varying highp vec2 vTextureCoord;
varying highp vec3 vFragPos;
varying highp vec3 vNormal;

// Shadow map related variables
#define NUM_SAMPLES 60
#define BLOCKER_SEARCH_NUM_SAMPLES NUM_SAMPLES
#define PCF_NUM_SAMPLES NUM_SAMPLES
#define NUM_RINGS 10

#define EPS 1e-3
#define PI 3.141592653589793
#define PI2 6.283185307179586

#define CAMERA_SIZE 400.0
#define Light_SIZE 50.0
#define RESOLUTION 2048.0

uniform sampler2D uShadowMap;

varying vec4 vPositionFromLight;

highp float rand_1to1(highp float x ) { 
  // -1 -1
  return fract(sin(x)*10000.0);
}

highp float rand_2to1(vec2 uv ) { 
  // 0 - 1
	const highp float a = 12.9898, b = 78.233, c = 43758.5453;
	highp float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );
	return fract(sin(sn) * c);
}

//  Input the sampled texture and return the depth
float unpack(vec4 rgbaDepth) {
    const vec4 bitShift = vec4(1.0, 1.0/255.0, 1.0/(255.0*255.0), 1.0/(255.0*255.0*255.0));
    return dot(rgbaDepth, bitShift);
}

vec2 poissonDisk[NUM_SAMPLES];

void poissonDiskSamples( const in vec2 randomSeed ) {

  float ANGLE_STEP = PI2 * float( NUM_RINGS ) / float( NUM_SAMPLES );
  float INV_NUM_SAMPLES = 1.0 / float( NUM_SAMPLES );

  float angle = rand_2to1( randomSeed ) * PI2;
  float radius = INV_NUM_SAMPLES;
  float radiusStep = radius;

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( cos( angle ), sin( angle ) ) * pow( radius, 0.75 );
    radius += radiusStep;
    angle += ANGLE_STEP;
  }
}

void uniformDiskSamples( const in vec2 randomSeed ) {

  float randNum = rand_2to1(randomSeed);
  float sampleX = rand_1to1( randNum ) ;
  float sampleY = rand_1to1( sampleX ) ;

  float angle = sampleX * PI2;
  float radius = sqrt(sampleY);

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( radius * cos(angle) , radius * sin(angle)  );

    sampleX = rand_1to1( sampleY ) ;
    sampleY = rand_1to1( sampleX ) ;

    angle = sampleX * PI2;
    radius = sqrt(sampleY);
  }
}

float getDynamicBias(float biasParaAdd){
  vec3 lightDir = normalize(uLightPos - vFragPos);
  vec3 normal = normalize(vNormal);
  float c = 0.005;
  c = biasParaAdd + c;
  float bias = max(c * (1.0 - dot(normal, lightDir)), c);
  return bias;
}

// ShadingPoint to Light. Try to find the shadowMap's range
// zReceiver 是接收光照表面点的深度值
float findBlocker( sampler2D shadowMap, vec2 uv, float zReceiver ) {
  // 声明参数
  float blockerTinyCount = 0.0;
  float totalBlockerDepth = 0.0;
  // 第一步 先对uv进行泊松圆盘采样
  poissonDiskSamples(uv);
  for(int i=0; i< NUM_SAMPLES; i++)
  {
    float blockerSearch = zReceiver * Light_SIZE / CAMERA_SIZE / 2.0;
    // poissonDisk[i] * 5.0/ RESOLUTION;
    vec4 shadowColor = texture2D(shadowMap, uv + poissonDisk[i] * blockerSearch);
    float shadowDepth = unpack(shadowColor);

    // 该微小部分在阴影中
    if(zReceiver - getDynamicBias(0.0) > shadowDepth){
      blockerTinyCount++;
      totalBlockerDepth += shadowDepth;
    }
  }
  // 不在范围内return 1;
  //if(blockerTinyCount - float(NUM_SAMPLES) <= EPS) return 1.0;
  if(blockerTinyCount == 0.0) return -1.0;
  // 最后返回平均深度
  float averageBlockerDepth = totalBlockerDepth / blockerTinyCount;
	return averageBlockerDepth;
}

float PCF(sampler2D shadowMap, vec4 coords, float filterSize) {
  // 先声明一下要用到的变量:最终返回的阴影值（可见度）visibility
  // 分辨率（此处在engine的定义为2048）
  // 采样所用的宽度filterSize:这里先采样5*5范围
  float visibility = 0.0;

  // 第一件事 先对这个坐标进行采样
  poissonDiskSamples(coords.xy);

  // 第二件事 对每个采样点周围区域进行深度值比较并且累加
  for(int i =0; i< PCF_NUM_SAMPLES; i++)
  {
    // 每个点乘上其覆盖一片区域内的uv范围
    vec2 uvOffset = poissonDisk[i] * filterSize/RESOLUTION;
    vec4 shadowColor = texture2D(shadowMap, coords.xy + uvOffset);
    float shadowDepth = unpack(shadowColor);
    // 计算当前点,同样加上bias
    float acc = coords.z < shadowDepth + getDynamicBias(0.0) ? 1.0:0.0;
    visibility += acc;
  }

  // 最后返回平均值
  return visibility / float(PCF_NUM_SAMPLES);
}

float PCSS(sampler2D shadowMap, vec4 coords){
  // 先进行泊松圆盘采样坐标的记录
  poissonDiskSamples(coords.xy);
  // STEP 1: avgblocker depth
  // 找到Blocker Distance
  float dBlocker = findBlocker(shadowMap, coords.xy, coords.z);
  //没有Blocker，直接可见
  if(dBlocker == -1.0) return 1.0;
  // STEP 2: penumbra size
  // 利用 STEP1 的BlockerDistance计算 filterSize
  float DReceiver = coords.z - dBlocker;
  float WPenumbra = DReceiver * Light_SIZE / dBlocker;
  // STEP 3: filtering
  // PCF传入 STEP2 计算的filterSize，5.0为魔法数字（可调整效果）
  return PCF(shadowMap, coords, WPenumbra);
}

// shadowCoord: xy is the texcoord, z records the cur depth and w is used for NDC division
float useShadowMap(sampler2D shadowMap, vec4 shadowCoord){
  vec4 shadowColor = texture2D(shadowMap, shadowCoord.xy);
  float shadowDepth = unpack(shadowColor);
  float curDepth = shadowCoord.z;

  float bias = getDynamicBias(0.01);
  if(shadowDepth + bias <= curDepth)
  {
    // 在阴影中
    return 0.0;
  }
  
  // 不在阴影中
  return 1.0;
}

vec3 blinnPhong() {
  vec3 color = texture2D(uSampler, vTextureCoord).rgb;
  color = pow(color, vec3(2.2));

  vec3 ambient = 0.05 * color;

  vec3 lightDir = normalize(uLightPos);
  vec3 normal = normalize(vNormal);
  float diff = max(dot(lightDir, normal), 0.0);
  vec3 light_atten_coff =
      uLightIntensity / pow(length(uLightPos - vFragPos), 2.0);
  vec3 diffuse = diff * light_atten_coff * color;

  vec3 viewDir = normalize(uCameraPos - vFragPos);
  vec3 halfDir = normalize((lightDir + viewDir));
  float spec = pow(max(dot(halfDir, normal), 0.0), 32.0);
  vec3 specular = uKs * light_atten_coff * spec;

  vec3 radiance = (ambient + diffuse + specular);
  vec3 phongColor = pow(radiance, vec3(1.0 / 2.2));
  return phongColor;
}

void main(void) {

  float visibility;
  // NDC 转换
  vec3 shadowCoord = vPositionFromLight.xyz / vPositionFromLight.w;

  // 将 [-1, 1] 映射到 [0, 1]
  shadowCoord = (shadowCoord.xyz + 1.0) * 0.5; 


  //visibility = useShadowMap(uShadowMap, vec4(shadowCoord, 1.0));
  //visibility = PCF(uShadowMap, vec4(shadowCoord, 1.0));
  visibility = PCSS(uShadowMap, vec4(shadowCoord, 1.0));

  vec3 phongColor = blinnPhong();

  gl_FragColor = vec4(phongColor * visibility, 1.0);
  //gl_FragColor = vec4(phongColor, 1.0);
}