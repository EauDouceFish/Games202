class DirectionalLight {

    constructor(lightIntensity, lightColor, lightPos, focalPoint, lightUp, hasShadowMap, gl) {
        this.mesh = Mesh.cube(setTransform(0, 0, 0, 0.2, 0.2, 0.2, 0));
        this.mat = new EmissiveMaterial(lightIntensity, lightColor);
        this.lightPos = lightPos;
        this.focalPoint = focalPoint;
        this.lightUp = lightUp

        this.hasShadowMap = hasShadowMap;
        this.fbo = new FBO(gl);
        if (!this.fbo) {
            console.log("无法设置帧缓冲区对象");
            return;
        }
    }

    // 光源MVP切换到投影空间，在PhongMaterial计算，LoadObj中创建
    CalcLightMVP(translate, scale) {
        let lightMVP = mat4.create();
        let modelMatrix = mat4.create();
        let viewMatrix = mat4.create();
        let projectionMatrix = mat4.create();
        
        // Model都需要声明最终输出在哪个矩阵上
        // Model transform
        mat4.translate(modelMatrix, modelMatrix,translate);
        mat4.scale(modelMatrix, modelMatrix,scale);

        // lookAt内置方法接收matrix, position, target, up
        // View transform
        mat4.lookAt(viewMatrix, this.lightPos, this.focalPoint, this.lightUp);

        // ortho使用正交投影，同样接收left,right,bottom,top,near,far
        // Projection transform
        const left = -200;
        const right = -left;
        const top = 200;
        const bottom = -top;
        const near = 0.01;
        const far = 1000;

        mat4.ortho(projectionMatrix,left,right,bottom,top,near,far);

        mat4.multiply(lightMVP, projectionMatrix, viewMatrix);
        mat4.multiply(lightMVP, lightMVP, modelMatrix);

        return lightMVP;
    }
}
