#ifndef PDFdef
#define PDFdef

#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#define PI 3.14159265358979323846
#define MAX(a,b) (((a)-(b)<=0)?(b):(a))
#define MIN(a,b) (((a)-(b)<=0)?(a):(b))

int MAX_length = 2000;
int MAX_slices = 1000;
double Interval = double(MAX_length)/MAX_slices;

void const_alter (int max_l, int max_s) {
    MAX_length = max_l;
    MAX_slices = max_s;
    Interval = double(MAX_length)/MAX_slices;
}

class PDF {
    public:
        double L_c;
        bool not_exist = false;
        std::vector<double> PDFarray;
        int length (void) {
            return this->PDFarray.size();
        }
        double Func (double x) {          // from array to func       this func is checked
            int sub = x / Interval;
            if (sub >= this->PDFarray.size()-1) 
                return 0;
            return PDFarray[sub] + ((PDFarray[sub+1]-PDFarray[sub])/Interval * (x - sub*Interval));     // linear interpolation
        }
        void show_y_array (void) {
            std::cout << "y ";
            for (int i = 0; i < this->PDFarray.size(); i++) {
                std::cout << this->PDFarray[i] << " ";
            }
            std::cout << "\n";
        }
        void show_x_array (void) {
            std::cout << "x ";
            for (int i = 0; i < this->PDFarray.size(); i++) {
                std::cout << i*Interval << " ";
            }
            std::cout << "\n";
        }
        void output(std::ostream &OutPut) { // output PDF to stream (file)
            OutPut << Interval << "\n";
            OutPut << this->L_c << "\n";
            for (int i = 0; i < this->PDFarray.size(); i++) {
                OutPut << this->PDFarray[i] << " \t";
            }
            OutPut << "\n";
        }
        PDF (std::istream &InPut) {     // get PDF from stream (file)
            std::string line;
            if (std::getline(InPut, line)) {
                Interval = std::stod(line);
            }
            if (std::getline(InPut, line)) {
                this->L_c = std::stod(line);
            }
            if (std::getline(InPut, line)) {
                std::istringstream line_stream(line);
                double num;
                while (line_stream >> num) {
                    this->PDFarray.push_back(num);
                }
            }
        }
        PDF (std::vector<double>& data, double L_c) {   // construct a PDF from data
            this->PDFarray = data;
            this->L_c = L_c;
        }
        PDF (char type, double diam) {           // construct a domain PDF
            if (type == 'd') {
                diam = diam;                  // most of the diams in this code are actually radius, because of my stupidity
                double max_d = diam*1.05 / Interval;
                double min_d = diam*0.95 / Interval;
                std::vector<double> data(max_d, 0);
                for (int i = min_d; i < data.size(); i++) {
                    data[i] = 1;
                }
                this->L_c = diam*1.05;
                this->PDFarray = data;
            }
            else 
                std::cout << "\nwrong type\n";
        }
        PDF (char type, double L_c, double L_p) { // construct a linker PDF
            if (type == 'l') {
                std::vector<double> data(L_c/Interval, 0);
                for (int i = 0; i < data.size(); i++) {
                    double r = i*Interval;
                    double a = (r / L_c);
                    a = a*a;
                    double b = (-9*L_c/8/L_p/(1-a));
                    data[i] = 1.0/ std::pow((1-a), 9.0/2) * std::exp(b);
                }
                this->L_c = L_c;
                this->PDFarray = data;
                this->normalize();
            }
            else
                std::cout << "\nwrong type\n";
        }
        PDF (void) {
            this->not_exist = true;
        }

        PDF (char type, PDF* f, PDF* g) {  // two ways of getting new PDF
            if (type == 'm') {  // multiply way
                if (f->not_exist || g->not_exist) { // judge if not exists
                    if (!g->not_exist) {        // only f not exists
                        this->L_c = g->L_c;
                        this->PDFarray = g->PDFarray;
                        return;
                    }
                    if (!f->not_exist) {        // only g not exists
                        this->L_c = f->L_c;
                        this->PDFarray = f->PDFarray;
                        return;
                    }
                    this->not_exist = true;
                    return;
                }
                this->L_c = MIN(f->L_c, g->L_c);
                std::vector<double> data( (this->L_c)/Interval , 0);
                for (int i = 0; i < data.size(); i++) {
                    data[i] = f->PDFarray[i] * g->PDFarray[i];
                }
                this->PDFarray = data;
            }
            if (type == 'c') {  // convolution way
                if (f->not_exist || g->not_exist) {
                    //std::cout << "error in PDFdef.cpp!\n";
                    this -> not_exist = true;
                    return;
                }

                std::vector<double> data( (f->L_c+g->L_c)/Interval , 0);
                this->L_c = f->L_c + g->L_c;

                if (f->L_c < g->L_c) {  // make sure f is the longer one
                    PDF* temp = f;
                    f = g;
                    g = temp;
                }
                // std::cout << "\n convoluting f: ";
                // f->show_y_array();
                // std::cout << "\n convoluting g: ";
                // g->show_y_array();
                
                /*
                5May2024 update:
                I have to change the way of convolution, because the previous way is too slow.
                The dimension to calculate each data point is lowered by one.
                And it's easy to proof that, for two radius functions f'(r), g'(r) defined on R^3, the convolution of them is:
                    (for simplicity, we use f(r): R|->R and g(r): R|->R to represent f'(r) and g'(r). So each variable here is a real number.)
                    (f*g)(R) = \frac{2\pi}{R} * \int_{0}^{r_{max}} r*f(r) * \int_{|R-r|}^{R+r} x*g(x) \,dx  \,dr        (1)
                there is no variable r in the second integral, so we can compute the integral of xg(x) first.
                so that the time compexity to calculate each data point is O(n) instead of O(n^2).
                But for the case R=0, the formula above turns out to be in a form of 0/0, so we have to calculate it separately with a formula:
                    (f*g)(0) = 2\pi * \int_{0}^{r_{max}} 2*r^2*f(r)*g(r) \,dr     (2)
                
                Now we do the convolution by the following steps:
                    1. compute the integral of xg(x)
                    2. compute the convolution of f(r) and g(r) for R=0 with formula (2)
                    3. compute the convolution of f(r) and g(r) for R>0 with formula (1)
                */

                // compute the integral of xg(x)
                double processed_g[g->PDFarray.size()];
                processed_g[0] = 0;
                for (int i = 1; i < g->PDFarray.size(); i++) {      // i*Interval = x
                    processed_g[i] = processed_g[i-1] + (g->PDFarray[i] * i);
                    // the above line should be "processed_g[i] = processed_g[i-1] + (g->PDFarray[i] * i*Interval) * Interval;"
                    // the first Interval is the x, the second Interval is the dx
                    // but we have another normalization step, so we can ignore all of the constants here.
                }
                
                // compute the convolution of f(r) and g(r) for R=0 with formula (2)
                data[0] = 0;
                for (int j = 1; j < g->PDFarray.size(); j++) {    // j*Interval = r
                    data[0] += j*j * f->PDFarray[j] * g->PDFarray[j];
                    // the above line should be "data[0] += j*j*Interval*Interval * f->PDFarray[j] * g->PDFarray[j] * Interval;"
                    // the first two Interval is the r^2, the thrid Interval is the dr
                }
                data[0] *= 2;   // the 2 in integral of formula (2). we ignore the 2pi here, and do the same to the formula (1)

                // compute the convolution of f(r) and g(r) for R>0 with formula (1)
                for (int i = 1; i < data.size(); i++) { // i*Interval = R
                    double sum = 0;
                    for (int j = 1; j < f->PDFarray.size(); j++) { // j*Interval = r
                        int integral_upper = (i+j > g->PDFarray.size()-1) ? g->PDFarray.size()-1 : i+j;     // the upper bound of the integral, g(x) = 0 where x > g->PDFarray.size()-1
                        int integral_lower = (std::abs(i-j) > g->PDFarray.size()-1) ? g->PDFarray.size()-1 : std::abs(i-j);
                        if (i > j) {
                            sum += j * f->PDFarray[j] * (processed_g[integral_upper] - processed_g[integral_lower]);
                        }
                        else if (i < j) {
                            sum += j * f->PDFarray[j] * (processed_g[integral_upper] + processed_g[integral_lower]);
                        }
                        else {
                            sum += j * f->PDFarray[j] * processed_g[integral_upper];
                        }
                        // the above line should be "sum += j*Interval * f->PDFarray[j] * (processed_g[i+j] - processed_g[std::abs(i-j)]) * Interval;"
                        // the first Interval is the r, the second Interval is the dr
                        // but we have another normalization step, so we can ignore all of the constants here.
                    }
                    sum /= (i);    // the 2pi/R in formula (1), we ignore the 2pi here
                    data[i] = sum;
                }
                this->PDFarray = data;
                // std::cout << "\n convoluted: ";
                // this->show_y_array();
                this->normalize();
            }
        }
        PDF (char type, PDF* f, PDF* g, double D) { // D refers to antigen density on the surface
            if (type == 's') {  // surface merge 2 sides binding restrictions
                PDF* long_PDF = (f->L_c > g->L_c) ? f : g;  // the longer one
                PDF* short_PDF = (f->L_c > g->L_c) ? f : g;  // the shorter one
                std::vector<double> data( short_PDF->L_c/Interval , 0);
                
                PDF* normalized_long = new PDF(long_PDF->PDFarray, long_PDF->L_c);
                normalized_long->TwoDnormalize();
                
                // compute PDF for R=0
                data[0] = 0;

                // compute PDF for R>0
                for (int i = 1; i < data.size(); i++) {
                    double sum = 0;
                    for (int j = 1; j < normalized_long->PDFarray.size(); j++) {
                        sum += j * normalized_long->PDFarray[j];
                    }
                    sum *= Interval * short_PDF->PDFarray[i] * 2 * PI;
                    data[i] = sum;
                }
                this->PDFarray = data;
                this->L_c = short_PDF->L_c;
                delete normalized_long;   // free the memory
                //this->normalize();
            }
        }
        void normalize (void) { // 3D normalization of the PDF
            double sum = 0;
            for (int i = 1; i < this->PDFarray.size()-1; i++) {
                sum += this->PDFarray[i] *i*i;
            }
            //sum += this->PDFarray[0] *4*PI*0*Interval*0*Interval/2;
            sum += this->PDFarray.back()/2.0 *(this->PDFarray.size()-1)*(this->PDFarray.size()-1);
            sum *= Interval*4*PI*Interval*Interval;
            for (int i = 0; i < this->PDFarray.size(); i++) {
                this->PDFarray[i] /= sum;
            }
        }
        void TwoDnormalize (void) { // this is a 2D normalization, which is used in the 2D convolution of cell surface simulation. 
            double sum = 0;
            for (int i = 1; i < this->PDFarray.size()-1; i++) {
                sum += this->PDFarray[i] *i;
            }
            sum += this->PDFarray.back()/2.0 *(this->PDFarray.size()-1);
            sum *= Interval*2*PI*Interval;
            for (int i = 0; i < this->PDFarray.size(); i++) {
                this->PDFarray[i] /= sum;
            }
        }


    private:
        double convolution_func (double phi, double r, double R, PDF* conv_func1, PDF* conv_func2) {    // this func is not in use, because of the new 1D convolution method
            double var = std::sqrt(R*R + r*r - 2*R*r*std::cos(phi));
            return conv_func1->Func(r) * conv_func2->Func(var) * r*r * std::sin(phi);
        }
        double DBLintegrate_for_conv(double x_L_bound, double x_U_bound, double y_L_bound, double y_U_bound, double min_step, double R, PDF* f, PDF* g) {   // this func is not in use, because of the new 1D convolution method
            int x_times = int( (x_U_bound-x_L_bound)/min_step +1 );
            int y_times = int( (y_U_bound-y_L_bound)/min_step +1 );
            double x_min_step = (x_U_bound-x_L_bound) / x_times;
            double y_min_step = (y_U_bound-y_L_bound) / y_times;
            double sum = 0;
            for (int i = 1; i < x_times; i++) {     // center
                for (int j = 1; j < y_times; j++) {
                    sum += convolution_func(x_L_bound+i*x_min_step, y_L_bound+j*y_min_step, R, f, g);
                }
            }
            sum += convolution_func(x_L_bound, y_L_bound, R, f, g)/3.0 + convolution_func(x_U_bound, y_U_bound, R, f, g)/3.0; // four corners
            sum += convolution_func(x_L_bound, y_U_bound, R, f, g)/6.0 + convolution_func(x_U_bound, y_L_bound, R, f, g)/6.0;
            double bound_sum = 0;
            for (int i = 1; i < x_times; i++) {     // sides
                bound_sum += convolution_func(x_L_bound+i*x_min_step, y_L_bound, R, f, g);
                bound_sum += convolution_func(x_L_bound+i*x_min_step, y_U_bound, R, f, g);
            }
            for (int i = 1; i < y_times; i++) {
                bound_sum += convolution_func(x_L_bound, y_L_bound+i*y_min_step, R, f, g);
                bound_sum += convolution_func(x_U_bound, y_L_bound+i*y_min_step, R, f, g);
            }
            bound_sum *= 0.5;
            sum += bound_sum;
            return sum*x_min_step*y_min_step;
        }
};


#endif