#pragma once

#include "Vec3D.h"
#include "Matrix.h"
#include "Constants.h"

class Projection {

    static const double  EtoPJ2000_0[3];
    static const double  EtoPJ2000_1[3];
    static const double  EtoPJ2000_2[3];
    static const double* const EtoP[3];


public:

    static Vec3D projectPosition(const Vec3D& target, const Vec3D& lookAt, const Vec3D& origin, const double fovAngle);
    static void project(Vec3D& position, Vec3D& velocity, Vec3D& lookAt, Vec3D& origin);
    static Vec3D projectVelocity(const Vec3D& target, const Vec3D& lookAt, const double fovAngle);


    //dimensions: [in/out] position: [pc] velocity: [km/s]
    static void GCAtoLSR(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with LSRtoGCA
    //dimensions: [in/out] position: [pc] velocity: [km/s]
    static void LSRtoGCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with GCAtoLSR

    static void GCAtoGCP(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut);
    static void GCPtoGCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut);

    //dimensions: [in/out] position: [pc] velocity: [km/s]
    static void LSRtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with HCAtoLSR
    static void HCAtoLSR(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with LSRtoHCA

    //epoch is J2000
    //dimensions: [in] position: [pc] velocity: [km/s] | [out] position: [pc,arcsec] velocity: [km/s,arcsec/yr]
    static void HCAtoHEQ(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with HEQtoHCA
    static void HEQtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with HCAtoHEQ
    static void HCAtoHGP(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with HGPtoHCA
    static void HGPtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut); //consistent with HCAtoHGP

    //static void HGPtoHEQ(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut);

};