
#ifndef CellSurfaceEffc
#define CellSurfaceEffc

#include "NodeDef.cpp"
#include "GraphStructure.cpp"
#include "PDFdef.cpp"
#include <vector>
    #include <fstream>
    #include <iostream>
#define aVOGADRO 0.0006022
// aVOGADRO = 6.022e23 / (1e9)^3,    (1e9)^3 is the constant when you transform Angstrom^3 to Liter

double SingleHemisphereVolumeRatio (double radius, double height) {
    // if the height is larger than the radius, the whole hemisphere is bound
    if (height >= radius) {
        return 1.0;
    }

    double Vcap = PI/3.0 * (height*height*height + 2.0*radius*radius*radius - 3.0*radius*radius*height);
    double Vhemisphere = 2.0/3.0 * PI * radius*radius*radius;
    return (Vhemisphere - Vcap) / Vhemisphere;
}

double RegularizedSpheicalWedgeVolume (double radius, double a, double c) {
    // derive b
    double b = sqrt(radius*radius - a*a - c*c);

    // get volume of a regularized spherical wedge
    double term1 = a*b*c / 3.0;
    double term2 = a * (a*a/3.0 - radius*radius) * std::atan(b/c);
    double term3 = 2.0/3.0 * radius*radius*radius * std::atan(b*a / (radius*c));
    return (term1 + term2 + term3);
}

double IntersectionBoundVolumeRatio (double distance, double radius1, double radius2, double height)
{
    // if the distance is larger than the sum of the radius, the two spheres are not intersecting
    if (distance >= radius1 + radius2) {
        return 0.0;
    }

    // if the distance is smaller than the difference of the radius, one sphere is completely inside the other. In this case, the volume ratio is the same as the single hemisphere volume ratio
    if (distance <= std::abs(radius1 - radius2)) {
        return SingleHemisphereVolumeRatio(std::min(radius1, radius2), height);
    }

    // get ceilling from helen formula
    double s = (distance + radius1 + radius2) / 2.0;
    double area = sqrt(s * (s - distance) * (s - radius1) * (s - radius2));
    double ceilling = 2.0 * area / distance;

    // if the ceilling is higher than the height, all of the volume is bound
    if (ceilling <= height) {
        return 1.0;
    }

    // derive a1 and a2 by solving triangles
    double a1 = sqrt(radius1*radius1 - ceilling*ceilling);
    double a2 = sqrt(radius2*radius2 - ceilling*ceilling);

    // get volume of intersections of two hemispheres
    double volume1 = PI/3.0 * (a1*a1*a1 + 2.0*radius1*radius1*radius1 - 3.0*radius1*radius1*a1);
    double volume2 = PI/3.0 * (a2*a2*a2 + 2.0*radius2*radius2*radius2 - 3.0*radius2*radius2*a2);
    double intersection_volume = (volume1 + volume2) / 2.0;   // divided by 2 because they are hemispheres

    // get volumes of general spheical wedges, which gives the total unbound volume
    double regularized_wedge1_1 = RegularizedSpheicalWedgeVolume(radius1, a1, height);
    double regularized_wedge1_2 = RegularizedSpheicalWedgeVolume(radius1, height, a1);
    double regularized_wedge2_1 = RegularizedSpheicalWedgeVolume(radius2, a2, height);
    double regularized_wedge2_2 = RegularizedSpheicalWedgeVolume(radius2, height, a2);
    double volume_unbound = regularized_wedge1_1 + regularized_wedge1_2 + regularized_wedge2_1 + regularized_wedge2_2;

    return (intersection_volume - volume_unbound) / intersection_volume;
}

double get_Pmeet (double radius1, double radius2, double height, double interval, double density) {
    double P_on_surface = 0.0;
    // single end bound
    if (radius2 <= 0) {
        P_on_surface =  SingleHemisphereVolumeRatio(radius1, height);
    }
    // double end bound
    else {
        for (double distance = interval; distance <= radius1 + radius2; distance += interval) {
            double ratio = IntersectionBoundVolumeRatio(distance, radius1, radius2, height);
            P_on_surface += ratio * distance;
        }
        P_on_surface /= (radius1 + radius2) * (radius1 + radius2) / 2.0;
        P_on_surface *= interval;   // times dr, the differential of distance
    }
    double P_antigen_on_surface = density;
    double P = P_on_surface * P_antigen_on_surface;

    return P;
}

int get_True_type(int n, int total_num_of_ligand, int receptor_length) { // get the ture type, so that we can get its parameters from Ligand_set
    if (n <= total_num_of_ligand) {
        return n;
    }
    else {
        return (n-1-total_num_of_ligand)/(receptor_length-1) +1;
    }
}

// effective concentration computing specifically on cell surface
double CellSurface_Effc (G_microstate* microstate1, G_microstate* microstate2, std::vector<G_ligand>& Ligand_set, std::vector<double> Density_surfaceLig) {
                        // Density_surfacelig is the array of the density of ligands on the surface of the cell (angstrom^-2)
    double effective_concentration = -2;    // -2 means the effective concentration is exceptionally uncomputed
    int binding_site = -1;
    int left_bind = -1;
    int right_bind = -1;
    
    for (int i = 0; i < microstate1->holder.size(); i++) {
        if (microstate1->holder[i].isVacant())
        if (!microstate2->holder[i].isVacant()) {   // i is the binding site of receptor in this transformation
            binding_site = i;
            break;
        }
    }

    // search for left and right binding site
    for (int i = binding_site-1; i >= 0; i--) {
        if (!microstate2->holder[i].isVacant()) {
            left_bind = i;
            break;
        }
    }
    for (int i = binding_site+1; i < microstate1->holder.size(); i++) {
        if (!microstate2->holder[i].isVacant()) {
            right_bind = i;
            break;
        }
    }

    // get contour length of left and right binding site
    double left_L_c = 0.0;
    double right_L_c = 0.0;
    if (left_bind == -1 && right_bind == -1) {  // didn't bind on left
            return -1;  // the first time to bind on cell surface, no need to compute effective concentration
    }
    if (right_bind != -1) {  // bound on right
        right_L_c = Ligand_set[0].domain[right_bind].diam * 1.05;
        for (int i = binding_site; i < right_bind; i++) {
            right_L_c += Ligand_set[0].linker[i].L_c;
            right_L_c += Ligand_set[0].domain[i].diam * 1.05;
        }
    }
    if (left_bind != -1) {  // bind on left
        left_L_c = Ligand_set[0].domain[left_bind].diam * 1.05;
        for (int i = left_bind; i < binding_site; i++) {
            left_L_c += Ligand_set[0].linker[i].L_c;
            left_L_c += Ligand_set[0].domain[i].diam * 1.05;
        }
    }

    int binding_antigen = get_True_type(microstate2->holder[binding_site].LigandType, Ligand_set.size()-1, microstate1->holder.size());

    double Pmeet = 0;
    if (left_bind == -1 || right_bind == -1) {  // only bind on one side
        Pmeet = get_Pmeet(((left_bind == -1) ? right_L_c : left_L_c), -1, 20*0.05, Interval, Density_surfaceLig[binding_antigen-1]);
    }
    else {  // bind on both sides
        Pmeet = get_Pmeet(left_L_c, right_L_c, 20*0.05, Interval, Density_surfaceLig[binding_antigen-1]);
    }

    return Pmeet /aVOGADRO;
}

#endif