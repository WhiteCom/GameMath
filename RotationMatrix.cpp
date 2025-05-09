#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define max(num1, num2) (num1) > (num2) ? (num1) : (num2)
#define min(num1, num2) (num1) < (num2) ? (num1) : (num2)

struct Vector3
{
    float x, y, z;

    Vector3(float _x = 0, float _y = 0, float _z = 0)
        : x(_x), y(_y), z(_z) {}
};

struct Quaternion
{
    float x, y, z, w;

    Quaternion(float _x = 0, float _y = 0, float _z = 0, float _w = 1)
        : x(_x), y(_y), z(_z), w(_w) {}
};

// 회전 행렬로부터 쿼터니언 추출
// 켄 슈메이크 알고리즘 사용 (회전행렬 -> 쿼터니언 추출)
Quaternion DecomposeQuaternionXYZWFromMatrix(float InMatrix[3][3])
{
    float root = 0.f;
    float trace = InMatrix[0][0] + InMatrix[1][1] + InMatrix[2][2];
    
    Quaternion result;

    float X = 0.f;
    float Y = 0.f;
    float Z = 0.f;
    float W = 0.f;

    if (trace > 0.f)
    {
        // W 요소를 구하고 나머지 X, Y, Z 계산
        root = sqrtf(trace + 1.f);
        W = 0.5f * root;
        root = 0.5f / root;
        X = (InMatrix[2][1] - InMatrix[1][2]) * root;
        Y = (InMatrix[0][2] - InMatrix[2][0]) * root;
        Z = (InMatrix[1][0] - InMatrix[0][1]) * root;
    }
    else
    {
        short i = 0;
        short j = 0;
        short k = 0;
        // X, Y, Z 중에서 가장 큰 요소 파악하기 (i에 저장)
        if (InMatrix[1][1] > InMatrix[0][0]) { i = 1; }
        if (InMatrix[2][2] > InMatrix[i][i]) { i = 2; }
        // i, j, k 순서 지정
        static const short next[3] = { 1, 2, 0 };
        j = next[i];
        k = next[j];
        // 가장 큰 요소의 값 구하기
        root = sqrtf(InMatrix[i][i] - InMatrix[j][j] - InMatrix[k][k] + 1.f);
        float* qt[3] = { &X, &Y, &Z };
        *qt[i] = 0.5f * root;

        root = 0.5f / root;
        // 나머지 두 요ㅗ소의 값 구하기
        *qt[j] = (InMatrix[i][j] + InMatrix[j][i]) * root;
        *qt[k] = (InMatrix[i][k] + InMatrix[k][i]) * root;
        // 마지막 W 값 구하기
        W = (InMatrix[j][k] - InMatrix[k][j]) * root;
    }

    result.x = X;
    result.y = Y;
    result.z = Z;
    result.w = W;
    
    return result;
}

Vector3 DecomposeEulerAngleXYZ(const float m[3][3])
{
    Vector3 result;

    // Y-Axis (Pitch)
    result.y = std::asin(-m[2][0]);
    
    // 짐벌락 체크 (cosTheta = 0)
    const float epsilon = 1e-6f;
    if (fabs(result.y - M_PI / 2) < epsilon) // Theta = 90 degree
    {
        result.z = 0.0f;
        result.x = std::atan2(m[0][1], m[0][2]);
    }
    else if (fabs(result.y + M_PI / 2) < epsilon) // Theta = -90 degree
    {
        result.x = 0.0f;
        result.z = std::atan2(-m[0][1], -m[0][2]);
    }
    else
    {
        // Yaw, Roll (Z-Axis, X-Axis)
        result.z = std::atan2(m[1][0], m[0][0]);
        result.x = std::atan2(m[2][1], m[2][2]);
    }

    // radian -> degree 변환
    result.x *= 180.0f / M_PI;
    result.y *= 180.0f / M_PI;
    result.z *= 180.0f / M_PI;

    return result;
}

// Quaternion -> Rotation Matrix -> Euler Angle
Vector3 QuaternionToEuler(const Quaternion& q)
{
    float m[3][3] = { {0.0f,}, };

    float x = q.x;
    float y = q.y;
    float z = q.z;
    float w = q.w;

    // Create Rotation Matrix
    m[0][0] = 1 - 2 * y * y - 2 * z * z;
    m[0][1] = 2 * x * y - 2 * w * z;
    m[0][2] = 2 * x * z + 2 * w * y;

    m[1][0] = 2 * x * y + 2 * w * z;
    m[1][1] = 1 - 2 * x * x - 2 * z * z;
    m[1][2] = 2 * y * z - 2 * w * x;

    m[2][0] = 2 * x * z - 2 * w * y;
    m[2][1] = 2 * y * z + 2 * w * x;
    m[2][2] = 1 - 2 * x * x - 2 * y * y;

    // Decompose Euler Angle From Rotation Matrix
    Vector3 result;
    result = DecomposeEulerAngleXYZ(m);

    return result;
}

Quaternion EulerToQuaternion(const Vector3& angle)
{
    float m[3][3] = { {0.0f,}, };

    // 각도(degre)를 라디안(radian)으로 변환
    float phi = angle.x * M_PI / 180.0f;
    float theta = angle.y * M_PI / 180.0f;
    float psi = angle.z * M_PI / 180.0f;

    // 아래는 시계방향으로 회전한다고 가정.

    // X-axis Rotation
    float cosPhi = std::cos(phi);
    float sinPhi = std::sin(phi);

    float Rx[3][3] = { {0.0f,}, };
    Rx[0][0] = 1.0f;
    Rx[0][1] = 0.0f;
    Rx[0][2] = 0.0f;

    Rx[1][0] = 0.0f;
    Rx[1][1] = cosPhi;
    Rx[1][2] = -sinPhi;

    Rx[2][0] = 0.0f;
    Rx[2][1] = sinPhi;
    Rx[2][2] = cosPhi;

    // Y-axis Rotation
    float cosTheta = std::cos(theta);
    float sinTheta = std::sin(theta);

    float Ry[3][3] = { {0.0f,}, };
    Ry[0][0] = cosTheta;
    Ry[0][1] = 0.0f;
    Ry[0][2] = -sinTheta;

    Ry[1][0] = 0.0f;
    Ry[1][1] = 1.0f;
    Ry[1][2] = 0.0f;

    Ry[2][0] = sinTheta;
    Ry[2][1] = 0.0f;
    Ry[2][2] = cosTheta;

    // Z-axis Rotation
    float cosPsi = std::cos(psi);
    float sinPsi = std::sin(psi);

    float Rz[3][3] = { {0.0f,}, };
    Rz[0][0] = cosPsi;
    Rz[0][1] = -sinPsi;
    Rz[0][2] = 0.0f;

    Rz[1][0] = sinPsi;
    Rz[1][1] = cosPsi;
    Rz[1][2] = 0.0f;

    Rz[2][0] = 0.0f;
    Rz[2][1] = 0.0f;
    Rz[2][2] = 1.0f;
    
    // 
    // case 1) X-Y-Z Matrix
    // m[0][0] = cosTheta*cosPsi 
    // m[0][1] = -cosTheta*sinPsi
    // m[0][2] = sinTheta
    // 
    // m[1][0] = cosPhi*sinPsi + sinPhi*sinTheta*cosPsi
    // m[1][1] = cosPhi*cosPsi - sinPhi*sinTheta*sinPsi
    // m[1][2] = -sinPhi*cosTheta
    // 
    // m[2][0] = sinPhi*sinPsi - cosPhi*sinTheta*cosPsi
    // m[2][1] = sinPhi*cosPsi + cosPhi*sinTheta*sinPsi
    // m[2][2] = cosPhi*cosTheta
    //
    
    //
    // case 1) X-Y-Z Matrix
    //

    // Rz * Ry
    float temp[3][3] = { {0.0f,}, };
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            temp[i][j] = 0.0f;
            for (int k = 0; k < 3; ++k)
            {
                temp[i][j] += Rz[i][k] * Ry[k][j];
            }
        }
    }
    
    // (Rz * Ry) * Rx
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            m[i][j] = 0.0f;
            for (int k = 0; k < 3; ++k)
            {
                m[i][j] += temp[i][k] * Rx[k][j];
            }
        }
    }
    
    Quaternion result;
    result = DecomposeQuaternionXYZWFromMatrix(m);
    return result;
}

int main()
{
    float m[12] = {
        1.000000, 0.000000, 0.000000,
        0.000000, 0.988713, 0.149822,
        0.000000, -0.149822, 0.988713,
        2.023960, 0.777257, -15.028656
    };

    // scale 추출
    float sx = sqrtf(m[0] * m[0] + m[3] * m[3] + m[6] * m[6]);
    float sy = sqrtf(m[1] * m[1] + m[4] * m[4] + m[7] * m[7]);
    float sz = sqrtf(m[2] * m[2] + m[5] * m[5] + m[8] * m[8]);

    printf("result Scale : [sx : %.6f], [sy : %.6f], [sz : %.6f]\n", sx, sy, sz);

    // 각 축 정규화
    float rotationMatrix[3][3] = { {0.0f, }, };

    // 회전 matrix 생성 (3 * 3)
    rotationMatrix[0][0] = m[0] / sx; rotationMatrix[0][1] = m[1] / sy; rotationMatrix[0][2] = m[2] / sz;
    rotationMatrix[1][0] = m[3] / sx; rotationMatrix[1][1] = m[4] / sy; rotationMatrix[1][2] = m[5] / sz;
    rotationMatrix[2][0] = m[6] / sx; rotationMatrix[2][1] = m[7] / sy; rotationMatrix[2][2] = m[8] / sz;

    // rotation 추출
    Quaternion q;
    float tr = m[0] + m[4] + m[8];

    Quaternion result = DecomposeQuaternionXYZWFromMatrix(rotationMatrix);
    printf("result Quaternion : [qx : %.6f], [qy : %.6f], [qz : %.6f], [qw : %.6f]\n", result.x, result.y, result.z, result.w);
    
    Vector3 eulerAngle_RotationMatrix = DecomposeEulerAngleXYZ(rotationMatrix);
    printf("result Euler Angle : [x : %.6f], [y : %.6f], [z : %.6f]\n", eulerAngle_RotationMatrix.x, eulerAngle_RotationMatrix.y, eulerAngle_RotationMatrix.z);

    // Quaternion -> Eulder
    Vector3 quaternionToEuler = QuaternionToEuler(result);
    printf("convert Quaternion to Euler : [x : %.6f], [y : %.6f], [z : %.6f]\n", quaternionToEuler.x, quaternionToEuler.y, quaternionToEuler.z);

    // Euler -> Quaternion
    Quaternion eulerToQuaternion = EulerToQuaternion(quaternionToEuler);
    printf("convert Euler To Quaternion : [x : %.6f], [y : %.6f], [z : %.6f], [w : %.6f]\n", eulerToQuaternion.x, eulerToQuaternion.y, eulerToQuaternion.z, eulerToQuaternion.w);

    return 0;
}