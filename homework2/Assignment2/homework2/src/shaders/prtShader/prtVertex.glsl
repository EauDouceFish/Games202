attribute vec3 aVertexPosition;
//attribute vec3 aNormalPosition;
// LightTransform矩阵保存在顶点
attribute mat3 aPrecomputeLT;

// 全局变量可以自动接收
uniform mat4 uModelMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uProjectionMatrix;
// 在材质中传递预计算的光照uniform
uniform mat3 uPrecomputeL[3];

// 传入片元着色器的参数
//varying highp mat3 vPrecomputeLT;
//varying highp vec3 vNormal;
varying highp vec3 vColor;

float LdotLT(mat3 PrecomputeL, mat3 PrecomputeLT) {
  vec3 L_0 = PrecomputeL[0];
  vec3 L_1 = PrecomputeL[1];
  vec3 L_2 = PrecomputeL[2];
  vec3 LT_0 = PrecomputeLT[0];
  vec3 LT_1 = PrecomputeLT[1];
  vec3 LT_2 = PrecomputeLT[2];
  return dot(L_0, LT_0) + dot(L_1, LT_1) + dot(L_2, LT_2);
}

void main(void) {
  // 无实际作用，避免aNormalPosition被优化后产生警告
  // vNormal = (uModelMatrix * vec4(aNormalPosition, 0.0)).xyz;

  for(int i = 0; i < 3; i++)
  {
    vColor[i] = LdotLT(aPrecomputeLT, uPrecomputeL[i]);
  }

  gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix * vec4(aVertexPosition, 1.0);
}