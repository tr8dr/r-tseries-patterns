/*
 * MIT License
 *
 * Copyright (c) 2015 Jonathan Shore
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

inline bool EQ (double a, double b) {
    return abs(a - b) < 1e7;
}

inline bool NE (double a, double b) {
    return abs(a - b) > 1e-7;
}

inline double max(double a, double b) {
    return a > b ? a : b;
}

inline double min(double a, double b) {
    return a < b ? a : b;
}

inline int max(int a, int b) {
    return a > b ? a : b;
}

inline int min(int a, int b) {
    return a < b ? a : b;
}

void apply_label (NumericVector& labels, int Istart, int Iend, double label) {
    for (int i = Istart ; i <= Iend ; i++)
        labels[i] = label;
}


// [[Rcpp::export]]
NumericVector label_direction (NumericVector y, double excursion, int Tinactive) {
    int len = y.size();
    NumericVector labels(len);

    if (len == 0)
        return labels;
    
    int Istart = 0;
    int Icursor = 0;

    int Imin = 0;
    int Imax = 0;
    double Vmin = y[0];
    double Vmax = y[0];
    double Vprior = y[0];

    while (Icursor < len) {
        double v = y[Icursor];
        Vprior = v;

        // determine whether there has been a retracement, requiring a split
        if ((Vmax - Vmin) >= excursion && Imin > Imax && (v - Vmin) >= excursion) {
            apply_label (labels, Istart, Imax-1, 0.0);
            apply_label (labels, Imax, Imin, -1.0);
            Istart = Imin;
            Imax = Icursor;
            Vmax = v;
        }
        else if ((Vmax - Vmin) >= excursion && Imax > Imin && (Vmax - v) >= excursion) {
            apply_label (labels, Istart, Imin-1, 0.0);
            apply_label (labels, Imin, Imax, +1.0);
            Istart = Imax;
            Imin = Icursor;
            Vmin = v;            
        }

        // check for "inactive" period where price has not progressed since latest min/max (upward direction)
        else if (Imax > Imin && (Icursor - Imax) >= Tinactive && (y[Icursor] - Vmax) < excursion) {
            if ((Vmax - Vmin) >= excursion) {
                apply_label (labels, Istart, Imin-1, 0.0);
                apply_label (labels, Imin, Imax, +1.0);
                apply_label (labels, Imax+1, Icursor, 0.0);
            } else
                apply_label (labels, Istart, Icursor, 0.0);
                
            Istart = Icursor;
            Imax = Icursor;
            Imin = Icursor;
            Vmax = v;
            Vmin = v;
        }

        // check for "inactive" period where price has not progressed since latest min/max (downward direction)
        else if (Imin > Imax && (Icursor - Imin) >= Tinactive && (Vmin - y[Icursor]) < excursion) {
            if ((Vmax - Vmin) >= excursion) {
                apply_label (labels, Istart, Imax-1, 0.0);
                apply_label (labels, Imax, Imin, -1.0);
                apply_label (labels, Imin+1, Icursor, 0.0);
            } else
                apply_label (labels, Istart, Icursor, 0.0);
                
            Istart = Icursor;
            Imax = Icursor;
            Imin = Icursor;
            Vmax = v;
            Vmin = v;
        }

        // adjust local max
        if (v >= Vmax)
            { Imax = Icursor;  Vmax = v; }

        // adjust local min
        if (v <= Vmin)
            { Imin = Icursor;  Vmin = v; }

        Icursor = Icursor + 1;
    }

    // finish end
    if ((Vmax - Vmin) >= excursion && Imin > Imax) {
        apply_label (labels, Istart, Imax-1, 0.0);
        apply_label (labels, Imax, Imin, -1.0);
        apply_label (labels, Imin+1, Icursor-1, 0.0);
    }
    else if ((Vmax - Vmin) >= excursion && Imax > Imin) {
        apply_label (labels, Istart, Imin-1, 0.0);
        apply_label (labels, Imin, Imax, +1.0);
        apply_label (labels, Imax+1, Icursor-1, 0.0);
    } else
        apply_label (labels, Istart, Icursor-1, 0.0);

    return labels;
}


// [[Rcpp::export]]
NumericVector filter_direction (NumericVector y, NumericVector vdir, double excursion) {
    int len = y.size();
    NumericVector ndir(len);

    if (len == 0)
        return ndir;

    int Ipos = 0;
    while (Ipos < len) {
        double dir = vdir[Ipos];
        if (dir == 0.0) {
            ndir[Ipos++] = 0.0;
            continue;
        }

        // locate end of region
        int Istart = Ipos;
        int Iend = Ipos;
        while (Iend < len && vdir[Iend] == dir) Iend++;
        Iend--;

        // maximum extent
        int Imaxfwd = Istart;
        int Imaxback = Iend;
        double Vmaxfwd = 0.0;
        double Vmaxback = 0.0;

        // determine ols in the forward direction
        double fExy = 0.0;
        double fExx = 0.0;
        double fEx = 0.0;
        double fEy = 0.0;

        for (int i = Istart ; i <= Iend ; i++) {
            double Xc = i - Istart;
            double Yc = y[i];
            fExy += Xc*Yc;
            fExx += Xc*Xc;
            fEx += Xc;
            fEy += Yc;

            double beta = (fExy - fEx*fEy/ (Xc+1.0)) / (fExx - fEx*fEx/ (Xc+1.0));
            double distance = dir * beta * Xc;
            if (distance > Vmaxfwd)
                { Vmaxfwd = distance; Imaxfwd = i; }
        }

        // determine ols in the backward direction
        double bExy = 0.0;
        double bExx = 0.0;
        double bEx = 0.0;
        double bEy = 0.0;

        for (int i = Iend ; i >= Istart ; i--) {
            double Xc = Iend - i;
            double Yc = y[i];
            bExy += Xc*Yc;
            bExx += Xc*Xc;
            bEx += Xc;
            bEy += Yc;

            double beta = (bExy - bEx*bEy/ (Xc+1.0)) / (bExx - bEx*bEx/ (Xc+1.0));
            double distance = -dir * beta * Xc;
            if (distance > Vmaxback)
                { Vmaxback = distance; Imaxback = i; }
        }

        Ipos = Iend+1;

        // label forward region if meets size requirement
        if (Vmaxfwd >= excursion) {
            for (int i = Istart; i <= Imaxfwd; i++)
                ndir[i] = dir;
            for (int i = Imaxfwd+1; i < Imaxback; i++)
                ndir[i] = 0.0;
        } else {
            for (int i = Istart; i <= Imaxback; i++)
                ndir[i] = 0.0;
        }

        // label backward region is meets size requirement
        if (Vmaxback >= excursion) {
            for (int i = Imaxback; i <= Iend; i++)
                ndir[i] = dir;
        } else {
            for (int i = max(Imaxback,Imaxfwd+1); i <= Iend; i++)
                ndir[i] = 0.0;
        }
        
    }

    return ndir;
}


