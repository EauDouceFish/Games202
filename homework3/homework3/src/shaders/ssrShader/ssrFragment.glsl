#ifdef GL_ES
precision highp float;
#endif

uniform vec3 uLightDir;
uniform vec3 uCameraPos;
uniform vec3 uLightRadiance;
uniform sampler2D uGDiffuse;
uniform sampler2D uGDepth;
uniform sampler2D uGNormalWorld;
uniform sampler2D uGShadow;
uniform sampler2D uGPosWorld;

varying mat4 vWorldToScreen;
varying highp vec4 vPosWorld;

#define M_PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307
#define INV_PI 0.31830988618
#define INV_TWO_PI 0.15915494309
#define EPSILON 0.0000001
#define RAY_MARCH_STEP 0.05
#define MAX_STEP_NUM 150

float Rand1(inout float p) {
  p = fract(p * .1031);
  p *= p + 33.33;
  p *= p + p;
  return fract(p);
}

vec2 Rand2(inout float p) {
  return vec2(Rand1(p), Rand1(p));
}

float InitRand(vec2 uv) {
	vec3 p3  = fract(vec3(uv.xyx) * .1031);
  p3 += dot(p3, p3.yzx + 33.33);
  return fract((p3.x + p3.y) * p3.z);
}

vec3 SampleHemisphereUniform(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = uv.x;
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(1.0 - z*z);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = INV_TWO_PI;
  return dir;
}

vec3 SampleHemisphereCos(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = sqrt(1.0 - uv.x);
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(uv.x);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = z * INV_PI;
  return dir;
}

void LocalBasis(vec3 n, out vec3 b1, out vec3 b2) {
  float sign_ = sign(n.z);
  if (n.z == 0.0) {
    sign_ = 1.0;
  }
  float a = -1.0 / (sign_ + n.z);
  float b = n.x * n.y * a;
  b1 = vec3(1.0 + sign_ * n.x * n.x * a, sign_ * b, -sign_ * n.x);
  b2 = vec3(b, sign_ + n.y * n.y * a, -n.y);
}

vec4 Project(vec4 a) {
  return a / a.w;
}

float GetDepth(vec3 posWorld) {
  float depth = (vWorldToScreen * vec4(posWorld, 1.0)).w;
  return depth;
}

/*
 * Transform point from world space to screen space([0, 1] x [0, 1])
 *
 */
vec2 GetScreenCoordinate(vec3 posWorld) {
  vec2 uv = Project(vWorldToScreen * vec4(posWorld, 1.0)).xy * 0.5 + 0.5;
  return uv;
}

float GetGBufferDepth(vec2 uv) {
  float depth = texture2D(uGDepth, uv).x;
  if (depth < 1e-2) {
    depth = 1000.0;
  }
  return depth;
}

vec3 GetGBufferNormalWorld(vec2 uv) {
  vec3 normal = texture2D(uGNormalWorld, uv).xyz;
  return normal;
}

vec3 GetGBufferPosWorld(vec2 uv) {
  vec3 posWorld = texture2D(uGPosWorld, uv).xyz;
  return posWorld;
}

float GetGBufferuShadow(vec2 uv) {
  float visibility = texture2D(uGShadow, uv).x;
  return visibility;
}

vec3 GetGBufferDiffuse(vec2 uv) {
  vec3 diffuse = texture2D(uGDiffuse, uv).xyz;
  diffuse = pow(diffuse, vec3(2.2));
  return diffuse;
}

bool RayMarch(vec3 ori, vec3 dir, out vec3 hitPos) {
  vec3 oneStep = normalize(dir) * RAY_MARCH_STEP;
  vec3 curPos = ori;
  // 当前存储的深度与当前对应的深度
  float uvGBufferDepth = 0.0;
  float curDepth = 0.0;
  for(int curStep = 0; curStep < MAX_STEP_NUM; curStep++)
  {
    // 每次步进一些距离
    curPos = curPos + oneStep;
    // 获得当前步进到位置的屏幕uv
    vec2 uvCurScreen = GetScreenCoordinate(curPos);
    // 用worldPos获取当前深度，用屏幕uv获取当前GBuffer深度
    curDepth = GetDepth(curPos);
    uvGBufferDepth = GetGBufferDepth(uvCurScreen);
    // 作比较
    if(curDepth - uvGBufferDepth > EPSILON)
    {
      hitPos = curPos;
      return true;
    }
  }
  return false;
}

/*
 * Evaluate diffuse bsdf value.
 *
 * wi, wo are all in world space.
 * uv is in screen space, [0, 1] x [0, 1].
 * 接收入射光，返回颜色
 */
vec3 EvalDiffuse (vec3 wi, vec2 uv) {
  vec3 L = vec3(0.0); 
  vec3 rho = GetGBufferDiffuse(uv);
  vec3 worldNormal = normalize(GetGBufferNormalWorld(uv));
  L = rho * INV_PI * max(0.0, dot(worldNormal, wi));
  return L;
}


/*
 * Evaluate directional light with shadow map
 * uv is in screen space, [0, 1] x [0, 1].
 *
 */
vec3 EvalDirectionalLight(vec2 uv) {
  vec3 Le = vec3(0.0);
  float visibility = GetGBufferuShadow(uv);
  Le = uLightRadiance * visibility;
  return Le;
}

/*
  Eval the reflection result using RayMarch
*/
vec3 EvalReflection(vec3 wi, vec2 uv)
{
  //vec3 ori = GetGBufferPosWorld(uv);
  vec3 ori = vPosWorld.xyz;
  vec3 worldNormal = GetGBufferNormalWorld(uv);
  vec3 dir = reflect(wi, worldNormal);
  vec3 hitPos = vec3(0.0);
  if(RayMarch(ori,dir, hitPos))
  {
    vec2 uvScreenHit = GetScreenCoordinate(hitPos);
    vec3 colorHit = GetGBufferDiffuse(uvScreenHit);
    return colorHit;
  }
  return vec3(0.0);
}


#define SAMPLE_NUM 1
void main() {
  float s = InitRand(gl_FragCoord.xy);
  vec3 wi = normalize(uLightDir);
  // wo为目标点到相机，为reflection的-wi
  vec3 wo = normalize(uCameraPos - vPosWorld.xyz);
  vec3 L = vec3(0.0);
  // 当前着色点所在的屏幕uv
  vec2 screen_uv = GetScreenCoordinate(vPosWorld.xyz);
  vec3 normal = GetGBufferNormalWorld(screen_uv);
  
  vec3 fr_and_cos = EvalDiffuse(wi, screen_uv);
  vec3 Li_and_view = EvalDirectionalLight(screen_uv);
  vec3 reflectionRayMarch = EvalReflection(-wo, screen_uv);

  L = fr_and_cos * Li_and_view;
  // 在此使用蒙特卡洛方法+SSR计算间接光照
  vec3 L_indirect = vec3(0.0);
  vec3 b1, b2;
  LocalBasis(normal, b1, b2);
  for(int i=0; i< SAMPLE_NUM; i++)
  {
    float pdf;
    vec3 direction = SampleHemisphereCos(s, pdf);
    direction = normalize(mat3(b1, b2, normal) * direction);
    vec3 hitPos;
    if(RayMarch(vPosWorld.xyz, direction, hitPos))
    {
      vec2 uvHitPosScreen = GetScreenCoordinate(hitPos);
      L_indirect += EvalDiffuse(direction, screen_uv)/pdf * EvalDiffuse(wi, uvHitPosScreen) * EvalDirectionalLight(uvHitPosScreen);
    }
  }
  L_indirect /= float(SAMPLE_NUM);

  //L += reflectionRayMarch;
  L = L + L_indirect;

  //初始结果
  //L = GetGBufferDiffuse(screen_uv);
  vec3 color = pow(clamp(L, vec3(0.0), vec3(1.0)), vec3(1.0 / 2.2));
  gl_FragColor = vec4(vec3(color.rgb), 1.0);
}
