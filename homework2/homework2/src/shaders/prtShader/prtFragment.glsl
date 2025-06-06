#ifdef GL_ES
precision mediump float;
#endif

varying highp vec3 vColor;

// toneMapping : toSRGB
vec3 toneMapping(vec3 color) {
    vec3 result;

    for (int i = 0; i < 3; ++i) {
        float value = color[i];

        if (value <= 0.0031308)
            result[i] = 12.92 * value;
        else
            result[i] = (1.0 + 0.055) * pow(value, 1.0 / 2.4) - 0.055;
    }

    return result;
}

void main(){
    vec3 finalColor = toneMapping(vColor);
    gl_FragColor = vec4(finalColor, 1.0);
}
