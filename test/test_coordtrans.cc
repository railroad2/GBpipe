#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

#define __MAIN__
#include "healpix_map.h"
#include "healpix_base.h"
#include "healpix_map_fitsio.h"
#include "rotmatrix.h"

using namespace std;

double deg2rad(double deg)
{
    return deg * M_PI/180.;
}

int test_compare()
{
    const int Nside = 128;
    const int Npix  = Nside * Nside * 12;

    Healpix_Base hb(Nside, RING, SET_NSIDE);
    Healpix_Map<int> m(Nside, RING, SET_NSIDE);
    Healpix_Map<int> m1(Nside, RING, SET_NSIDE);

    rotmatrix R1;
    R1.Make_CPAC_Euler_Matrix(deg2rad(60), deg2rad(30), deg2rad(0));
    cout << "Rotation matrix = \n" << R1 << endl;
    R1.Transpose();

    for (int i=0; i<Npix; i++) {
        m[i] = i;
    } 
    for (int i=0; i<Npix; i++) {
        vec3 vm1 = hb.pix2vec(i);
        vec3 vm  = R1.Transform(vm1);
        int idx = hb.vec2pix(vm);
        m1[i] = m[idx];
    }

    write_Healpix_map_to_fits("m_original.fits", m,  PLANCK_INT64);
    write_Healpix_map_to_fits("m_rotation.fits", m1, PLANCK_INT64);

    return 0;
}

int test_EG()
{
    const int Nside = 128;
    const int Npix  = Nside * Nside * 12;

    Healpix_Base hb(Nside, RING, SET_NSIDE);
    Healpix_Map<int> m(Nside, RING, SET_NSIDE);
    Healpix_Map<int> m1(Nside, RING, SET_NSIDE);

    rotmatrix R1(-0.054876, -0.873437, -0.483835, \
                  0.494109, -0.444830,  0.746982, \
                 -0.867666, -0.198076,  0.455984);
    R1.Transpose();

    for (int i=0; i<Npix; i++) {
        m[i] = i;
    } 
    for (int i=0; i<Npix; i++) {
        vec3 vm1 = hb.pix2vec(i);
        vec3 vm  = R1.Transform(vm1);
        int idx = hb.vec2pix(vm);
        m1[i] = m[idx];
    }

    write_Healpix_map_to_fits("m_original.fits", m,  PLANCK_INT64);
    write_Healpix_map_to_fits("m_rotation.fits", m1, PLANCK_INT64);

    return 0;

}

int test_rotmatrices()
{
    rotmatrix R1, R2, R3;
    R1.Make_CPAC_Euler_Matrix(deg2rad(0), deg2rad(30), deg2rad(30));
    R2.Make_CPAC_Euler_Matrix(deg2rad(0), deg2rad(30), deg2rad(0));
    R3.Make_CPAC_Euler_Matrix(deg2rad(0), deg2rad(0), deg2rad(30));

    cout << R1 << endl;
    cout << R2 << endl;
    cout << R3 << endl;
    
    return 0;
}

int main()
{
    //test_compare();
    //test_EG();
    test_rotmatrices();

    return 0;
}
