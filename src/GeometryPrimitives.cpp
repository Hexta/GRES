/******************************************************************************
 * Copyright (c) 2009-2016 Artur Molchanov <artur.molchanov@gmail.com>        *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 ******************************************************************************/

#define GL_GLEXT_PROTOTYPES
#include "GeometryPrimitives.h"

#include "Atoms.h"
#include "Coords3D.h"

#ifdef _WIN32
#include <windows.h>
#include <wingdi.h>
#endif

#include <GL/glu.h>
#include <GL/glext.h>

#include <cmath>

#define CACHE_SIZE 240
#define PI 3.1415926535897932385
#define SIN sinf
#define COS cosf
#define SQRT sqrtf

#ifdef _WIN32
#define glBindBufferARB pglBindBufferARB
#define glDeleteBuffersARB pglDeleteBuffersARB
#define glBufferDataARB pglBufferDataARB
#define glBufferSubDataARB pglBufferSubDataARB
#endif

namespace GeometryPrimitives {
void createSphere(double radius, int slices, int stacks, int& vSize1, int& vSize2, int& vSize3) {
#ifdef _WIN32
    PFNGLBINDBUFFERARBPROC pglBindBufferARB =
        reinterpret_cast<PFNGLBINDBUFFERARBPROC>(wglGetProcAddress("glBindBufferARB"));
    PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB =
        reinterpret_cast<PFNGLDELETEBUFFERSARBPROC>(wglGetProcAddress("glDeleteBuffersARB"));
    PFNGLBUFFERDATAARBPROC pglBufferDataARB =
        reinterpret_cast<PFNGLBUFFERDATAARBPROC>(wglGetProcAddress("glBufferDataARB"));
    PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB =
        reinterpret_cast<PFNGLBUFFERSUBDATAARBPROC>(wglGetProcAddress("glBufferSubDataARB"));
#endif
    const auto SIZE_OF_FLOAT = sizeof(float);

    GLint i, j;
    GLfloat sinCache1a[CACHE_SIZE];
    GLfloat cosCache1a[CACHE_SIZE];
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
    GLboolean needCache3;
    GLint start, finish;

    if (slices >= CACHE_SIZE) slices = CACHE_SIZE - 2;
    if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE - 1;
    if (slices < 2 || stacks < 1 || radius < 0.0) {
        return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */

    needCache3 = GL_TRUE;

    for (i = 0; i < slices; ++i) {
        angle = static_cast<GLfloat>(2 * PI * i / slices);
        sinCache1a[i] = static_cast<GLfloat>(SIN(angle));
        cosCache1a[i] = static_cast<GLfloat>(COS(angle));
    }

    for (j = 0; j <= stacks; ++j) {
        angle = static_cast<GLfloat>(PI * j / stacks);
        sinCache1b[j] = static_cast<GLfloat>(radius * SIN(angle));
        cosCache1b[j] = static_cast<GLfloat>(radius * COS(angle));
    }
    /* Make sure it comes to a point */
    sinCache1b[0] = 0;
    sinCache1b[stacks] = 0;

    if (needCache3) {
        for (i = 0; i < slices; ++i) {
            angle = static_cast<GLfloat>(2 * PI * (i - 0.5) / slices);
            sinCache3a[i] = static_cast<GLfloat>(SIN(angle));
            cosCache3a[i] = static_cast<GLfloat>(COS(angle));
        }
        for (j = 0; j <= stacks; ++j) {
            angle = static_cast<GLfloat>(PI * (j - 0.5) / stacks);
            sinCache3b[j] = static_cast<GLfloat>(SIN(angle));
            cosCache3b[j] = static_cast<GLfloat>(COS(angle));
        }
    }

    sinCache1a[slices] = sinCache1a[0];
    cosCache1a[slices] = cosCache1a[0];

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

    // GL_TRIANGLE_FAN

    Coords3DList norm1;
    Coords3DList vertex1;

    vertex1.reserve(slices + 1);
    norm1.reserve(slices + 1);

    Coords3D v1 = {0.0, 0.0, static_cast<float>(radius)};
    vertex1.push_back(v1);
    v1.normalize();
    norm1.push_back(v1);
    for (i = slices; i >= 0; --i) {
        Coords3D n1 = {sinCache3a[i] * sintemp3,
                       cosCache3a[i] * sintemp3, costemp3};
        norm1.push_back(n1);
        Coords3D v1 = {sintemp2* sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
        vertex1.push_back(v1);
    }

    /* High end next (j == stacks-1 iteration) */
    sintemp2 = sinCache1b[stacks - 1];
    zHigh = cosCache1b[stacks - 1];
    sintemp3 = sinCache3b[stacks];
    costemp3 = cosCache3b[stacks];
    // GL_TRIANGLE_FAN

    Coords3DList norm2;
    Coords3DList vertex2;

    vertex2.reserve(slices + 1);
    norm2.reserve(slices + 1);

    Coords3D v2 = {0.0, 0.0, static_cast<float>(-radius)};
    vertex2.push_back(v2);
    v2.normalize();
    norm2.push_back(v2);
    for (i = 0; i <= slices; ++i) {
        Coords3D n2 = {sinCache3a[i] * sintemp3, cosCache3a[i] * sintemp3, costemp3};
        norm2.push_back(n2);
        Coords3D v2 = {sintemp2* sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
        vertex2.push_back(v2);
    }

    Coords3DList norm3;
    Coords3DList vertex3;

    vertex3.reserve(2 * (finish - start) * slices);
    norm3.reserve(2 * (finish - start) * slices);

    for (j = start; j < finish; ++j) {
        zLow = cosCache1b[j];
        zHigh = cosCache1b[j + 1];
        sintemp1 = sinCache1b[j];
        sintemp2 = sinCache1b[j + 1];
        sintemp4 = sinCache3b[j + 1];
        costemp4 = cosCache3b[j + 1];

        // GL_QUAD_STRIP
        for (i = 0; i <= slices; ++i) {
            Coords3D v3 = {sintemp2* sinCache1a[i], sintemp2 * cosCache1a[i], zHigh};
            vertex3.push_back(v3);

            Coords3D n3 = {sinCache3a[i] * sintemp4, cosCache3a[i] * sintemp4, costemp4};
            norm3.push_back(n3);
            norm3.push_back(n3);

            Coords3D v3_ = {sintemp1* sinCache1a[i], sintemp1 * cosCache1a[i], zLow};
            vertex3.push_back(v3_);
        }
    }
    vSize1 = static_cast<int>(vertex1.size());

    glBindBufferARB(GL_ARRAY_BUFFER, 1);
    glBufferDataARB(GL_ARRAY_BUFFER, (vSize1 + norm1.size()) * 3 * SIZE_OF_FLOAT, 0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vSize1 * 3 * SIZE_OF_FLOAT, &vertex1[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vSize1 * 3 * SIZE_OF_FLOAT,
        norm1.size() * 3 * SIZE_OF_FLOAT, &norm1[0].x);

    vSize2 = static_cast<int>(vertex2.size());

    glBindBufferARB(GL_ARRAY_BUFFER, 2);
    glBufferDataARB(GL_ARRAY_BUFFER, (vSize2 + norm2.size()) * 3 * SIZE_OF_FLOAT,
        0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vSize2 * 3 * SIZE_OF_FLOAT,
        &vertex2[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vSize2 * 3 * SIZE_OF_FLOAT,
        norm2.size() * 3 * SIZE_OF_FLOAT, &norm2[0].x);

    vSize3 = static_cast<int>(vertex3.size());

    glBindBufferARB(GL_ARRAY_BUFFER, 3);
    glBufferDataARB(GL_ARRAY_BUFFER, (vSize3 + norm3.size()) * 3 * SIZE_OF_FLOAT,
        0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, vSize3 * 3 * SIZE_OF_FLOAT,
        &vertex3[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, vSize3 * 3 * SIZE_OF_FLOAT,
        norm3.size() * 3 * SIZE_OF_FLOAT, &norm3[0].x);

    glBindBufferARB(GL_ARRAY_BUFFER, 0);
    vertex1.clear();
    norm1.clear();
    vertex2.clear();
    norm2.clear();
    vertex3.clear();
    norm3.clear();
}
} // namespace GeometryPrimitives
