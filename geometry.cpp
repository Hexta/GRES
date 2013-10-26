#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include "geometry.h"

#define CACHE_SIZE 240
#define PI 3.1415926535897932385
#define SIN sinf
#define COS cosf
#define SQRT sqrtf

#ifdef _WIN32
#include <windows.h>
#define glBindBufferARB           pglBindBufferARB
#define glDeleteBuffersARB        pglDeleteBuffersARB
#define glBufferDataARB           pglBufferDataARB
#define glBufferSubDataARB        pglBufferSubDataARB
#endif

void
createAtomsAndBondes(surface3D &surface, vector<atomType> &surfAtoms,
                     atomsCoords &cellAts, float xs, float ys, float zs,
                     int z_min, float scaling, vector<atomName> &atNames_,
                     Bonds &outBonds) {

    int name = 0;
    for (auto &surfAtom : surfAtoms) {

        const int x = surfAtom.x;
        const int y = surfAtom.y;
        const int z = surfAtom.z;
        const unsigned char a = surfAtom.type;

        if (!surface[z][y][x][a].neighbours.empty()) {
            float x0 = scaling * x*xs;
            float y0 = scaling * y*ys;
            float z0 = -scaling * (z - z_min) * zs;

            float xA = x0 + scaling * cellAts[a].x;
            float yA = y0 + scaling * cellAts[a].y;
            float zA = z0 - scaling * cellAts[a].z;

            ++name;
            atomName temp = {name, x, y, z, xA, yA, zA, a, (short) surface[z][y][x][a].neighbours.size()};
            atNames_.push_back(temp);

            for (auto &nb : surface[z][y][x][a].neighbours) {
                float xNb = scaling * (xs * nb.x + cellAts[nb.type].x);
                float yNb = scaling * (ys * nb.y + cellAts[nb.type].y);
                float zNb = scaling * (-zs * nb.z - cellAts[nb.type].z);

                Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                outBonds.push_back(bondT);
            }
        }
    }

    for (auto &surfAtom : surfAtoms) {
        short x_ = surfAtom.x;
        short y_ = surfAtom.y;
        unsigned short z_ = surfAtom.z;
        unsigned char a_ = surfAtom.type;

        for (auto &nb : surface[z_][y_][x_][a_].neighbours) {
            //unsigned short a = (*surfAtoms)[i].type;

            short x = nb.x;
            short y = nb.y;
            short z = nb.z;
            unsigned char a = nb.type;

            if (x > 1 && x < surface[z_][y_].size() - 2 && y > 1 && y < surface[z_].size() - 2)
                if (!surface[z][y][x][a].deleted) {

                    float x0 = scaling * x*xs;
                    float y0 = scaling * y*ys;
                    float z0 = -scaling * (z - z_min) * zs;

                    float xA = x0 + scaling * cellAts[a].x;
                    float yA = y0 + scaling * cellAts[a].y;
                    float zA = z0 - scaling * cellAts[a].z;

                    ++name;
                    atomName temp = {name, x, y, z, xA, yA, zA, a,
                                     (short) surface[z][y][x][a].neighbours.size()};
                    atNames_.push_back(temp);

                    for (auto &nb_int : surface[z][y][x][a].neighbours) {
                        float xNb = scaling * (xs * nb_int.x + cellAts[nb_int.type].x);
                        float yNb = scaling * (ys * nb_int.y + cellAts[nb_int.type].y);
                        float zNb = scaling * (-zs * nb_int.z - cellAts[nb_int.type].z);

                        Bond bondT = {xA, yA, zA,
                                      xNb, yNb, zNb};
                        outBonds . push_back(bondT);
                    }
                }
        }
    }
}

void
normalize(float v[3]) {
    const GLfloat len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
    glNormal3fv(&v[0]);
}

void
normalize(float v[3], coords3D &out) {
    const GLfloat len = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
    out.x = v[0] / len;
    out.y = v[1] / len;
    out.z = v[2] / len;
}

coords3D
normalize(coords3D in) {
    GLfloat len = sqrt(pow(in.x, 2) + pow(in.y, 2) + pow(in.z, 2));
    return
    {
        in.x / len, in.y / len, in.z / len
    };
}

void
createSphere(GLdouble radius, GLint slices, GLint stacks, int* vSize1, int* vSize2, int* vSize3) {
#ifdef _WIN32
    PFNGLBINDBUFFERARBPROC pglBindBufferARB = (PFNGLBINDBUFFERARBPROC) wglGetProcAddress("glBindBufferARB");
    PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) wglGetProcAddress("glDeleteBuffersARB");
    PFNGLBUFFERDATAARBPROC pglBufferDataARB = (PFNGLBUFFERDATAARBPROC) wglGetProcAddress("glBufferDataARB");
    PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB = (PFNGLBUFFERSUBDATAARBPROC) wglGetProcAddress("glBufferSubDataARB");
#endif
    GLint i, j;
    GLfloat sinCache1a[CACHE_SIZE];
    GLfloat cosCache1a[CACHE_SIZE];
    GLfloat sinCache2a[CACHE_SIZE];
    GLfloat cosCache2a[CACHE_SIZE];
    GLfloat sinCache3a[CACHE_SIZE];
    GLfloat cosCache3a[CACHE_SIZE];
    GLfloat sinCache1b[CACHE_SIZE];
    GLfloat cosCache1b[CACHE_SIZE];
    GLfloat sinCache3b[CACHE_SIZE];
    GLfloat cosCache3b[CACHE_SIZE];
    GLfloat angle;
    GLfloat zLow, zHigh;
    GLfloat sintemp1 = 0.0, sintemp2 = 0.0, sintemp3 = 0.0, sintemp4 = 0.0;
    GLfloat costemp3 = 0.0, costemp4 = 0.0;
    GLboolean needCache2, needCache3;
    GLint start, finish;

    atomsCoords norm1, norm2, norm3;
    atomsCoords vertex1, vertex2, vertex3;

    if (slices >= CACHE_SIZE) slices = CACHE_SIZE - 1;
    if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE - 1;
    if (slices < 2 || stacks < 1 || radius < 0.0) {
        return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */
    needCache2 = needCache3 = GL_FALSE;

    needCache3 = GL_TRUE;

    for (i = 0; i < slices; ++i) {
        angle = 2 * PI * i / slices;
        sinCache1a[i] = SIN(angle);
        cosCache1a[i] = COS(angle);
        if (needCache2) {
            sinCache2a[i] = sinCache1a[i];
            cosCache2a[i] = cosCache1a[i];
        }
    }

    for (j = 0; j <= stacks; ++j) {
        angle = PI * j / stacks;
        sinCache1b[j] = radius * SIN(angle);
        cosCache1b[j] = radius * COS(angle);
    }
    /* Make sure it comes to a point */
    sinCache1b[0] = 0;
    sinCache1b[stacks] = 0;

    if (needCache3) {
        for (i = 0; i < slices; ++i) {
            angle = 2 * PI * (i - 0.5) / slices;
            sinCache3a[i] = SIN(angle);
            cosCache3a[i] = COS(angle);
        }
        for (j = 0; j <= stacks; ++j) {
            angle = PI * (j - 0.5) / stacks;
            sinCache3b[j] = SIN(angle);
            cosCache3b[j] = COS(angle);
        }
    }

    sinCache1a[slices] = sinCache1a[0];
    cosCache1a[slices] = cosCache1a[0];
    if (needCache2) {
        sinCache2a[slices] = sinCache2a[0];
        cosCache2a[slices] = cosCache2a[0];
    }
    if (needCache3) {
        sinCache3a[slices] = sinCache3a[0];
        cosCache3a[slices] = cosCache3a[0];
    }

    /* Do ends of sphere as TRIANGLE_FAN's (if not texturing)
     ** We don't do it when texturing because we need to respecify the
     ** texture coordinates of the apex for every adjacent vertex (because
     ** it isn't a constant for that point)
     */
    start = 1;
    finish = stacks - 1;

    /* Low end first (j == 0 iteration) */
    sintemp2 = sinCache1b[1];
    zHigh = cosCache1b[1];
    sintemp3 = sinCache3b[1];
    costemp3 = cosCache3b[1];

    //GL_TRIANGLE_FAN

    coords3D v1 = {0.0, 0.0, (float) radius};
    vertex1.push_back(v1);
    coords3D nT;
    normalize(&v1.x, nT);
    norm1.push_back(nT);
    for (i = slices; i >= 0; --i) {
        /*if (i != slices)*/
        {
            coords3D n1 = {sinCache3a[i + 1] * sintemp3,
                           cosCache3a[i + 1] * sintemp3, costemp3};
            norm1.push_back(n1);
        }
        coords3D v1 = {sintemp2 * sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
        vertex1.push_back(v1);
    }

    /* High end next (j == stacks-1 iteration) */
    sintemp2 = sinCache1b[stacks - 1];
    zHigh = cosCache1b[stacks - 1];
    sintemp3 = sinCache3b[stacks];
    costemp3 = cosCache3b[stacks];
    //GL_TRIANGLE_FAN

    coords3D v2 = {0.0, 0.0, (float) -radius};
    vertex2.push_back(v2);
    normalize(&v2.x, nT);
    norm2.push_back(nT);
    for (i = 0; i <= slices; ++i) {
        coords3D n2 = {sinCache3a[i] * sintemp3, cosCache3a[i] * sintemp3, costemp3};
        norm2.push_back(n2);
        coords3D v2 = {sintemp2 * sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
        vertex2.push_back(v2);
    }

    for (j = start; j < finish; ++j) {
        zLow = cosCache1b[j];
        zHigh = cosCache1b[j + 1];
        sintemp1 = sinCache1b[j];
        sintemp2 = sinCache1b[j + 1];
        sintemp4 = sinCache3b[j + 1];
        costemp4 = cosCache3b[j + 1];

        //GL_QUAD_STRIP
        for (i = 0; i <= slices; ++i) {
            coords3D v3 = {sintemp2 * sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
            vertex3.push_back(v3);

            coords3D n3 = {sinCache3a[i] * sintemp4, cosCache3a[i] * sintemp4, costemp4};
            norm3.push_back(n3);
            norm3.push_back(n3);

            coords3D v3_ = {sintemp1 * sinCache1a[i], sintemp1 * cosCache1a[i], zLow};
            vertex3.push_back(v3_);
        }
    }
    *vSize1 = vertex1.size();

    glBindBufferARB(GL_ARRAY_BUFFER, 1);
    glBufferDataARB(GL_ARRAY_BUFFER, (vertex1.size() + norm1.size())*3 * sizeof (float), 0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vertex1.size()*3 * sizeof (float), &vertex1[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vertex1.size()*3 * sizeof (float), norm1.size()*3 * sizeof (float), &norm1[0].x);

    *vSize2 = vertex2.size();

    glBindBufferARB(GL_ARRAY_BUFFER, 2);
    glBufferDataARB(GL_ARRAY_BUFFER, (vertex2.size() + norm2.size())*3 * sizeof (float),
                    0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vertex2.size()*3 * sizeof (float),
                       &vertex2[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vertex2.size()*3 * sizeof (float),
                       norm2.size()*3 * sizeof (float), &norm2[0].x);

    *vSize3 = vertex3.size();

    glBindBufferARB(GL_ARRAY_BUFFER, 3);
    glBufferDataARB(GL_ARRAY_BUFFER, (vertex3.size() + norm3.size())*3 * sizeof (float),
                    0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vertex3.size()*3 * sizeof (float),
                       &vertex3[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vertex3.size()*3 * sizeof (float),
                       norm3.size()*3 * sizeof (float), &norm3[0].x);

    glBindBufferARB(GL_ARRAY_BUFFER, 0);
    vertex1.clear();
    norm1.clear();
    vertex2.clear();
    norm2.clear();
    vertex3.clear();
    norm3.clear();
}

void
norm(coords3D &in) {
    GLfloat len = sqrt(pow(in.x, 2) + pow(in.y, 2) + pow(in.z, 2));
    coords3D temp = {in.x / len, in.y / len, in.z / len};
    in = temp;
}

coords3D
normcrossprod(coords3D v1, coords3D v2) {
    coords3D out;
    out.x = v1.y * v2.z - v1.z * v2.y;
    out.y = v1.z * v2.x - v1.x * v2.z;
    out.z = v1.x * v2.y - v1.y * v2.x;
    norm(out);
    return out;
}
