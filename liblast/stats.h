// stats.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010

/*
    This file is part of LibLast (liblast).

    LibLast is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LibLast is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LibLast.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef STATS_H
#define STATS_H

#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>

using namespace std;

template<int N, class T>
T nthPower(T x) {
    T ret = x;
    for (int i=1; i < N; ++i) {
        ret *= x;
    }
    return ret;
}

template<class T, int N>
struct SumDiffNthPower {
    SumDiffNthPower(T x) : mean_(x) { };
    T operator( )(T sum, T current) {
        return sum + nthPower<N>(current - mean_);
    }
    T mean_;
};

template<class T, int N, class Iter_T>
T nthMoment(Iter_T first, Iter_T last, T mean)  {
    size_t cnt = distance(first, last);
    return accumulate(first, last, T( ), SumDiffNthPower<T, N>(mean)) / cnt;
}

template<class T, class Iter_T>
T computeVariance(Iter_T first, Iter_T last, T mean) {
    return nthMoment<T, 2>(first, last, mean);
}

template<class T, class Iter_T>
T computeStdDev(Iter_T first, Iter_T last, T mean) {
    return sqrt(computeVariance(first, last, mean));
}

template<class T, class Iter_T>
T computeSkew(Iter_T begin, Iter_T end, T mean) {
    T m3 = nthMoment<T, 3>(begin, end, mean);
    T m2 = nthMoment<T, 2>(begin, end, mean);
    return m3 / (m2 * sqrt(m2));
}

template<class T, class Iter_T>
T computeKurtosisExcess(Iter_T begin, Iter_T end, T mean) {
    T m4 = nthMoment<T, 4>(begin, end, mean);
    T m2 = nthMoment<T, 2>(begin, end, mean);
    return m4 / (m2 * m2) - 3;
}

// AM
template<class T, class Iter_T>
T computeMedian(Iter_T begin, Iter_T end, T sum) {
    Iter_T s_it;
    T s_sum = 0;
    T median = 0;
    s_it = begin;
    while (true) {
        s_sum += (*s_it);
        if (abs(s_sum) >= abs(sum/2.0)) break;
        if (s_it == end) exit(1); // this shouldn't happen
        s_it++;
    }
    median = (*s_it);
    if (s_it!=begin) median = (median + (*(s_it--)))/2.0;
    return median;
}


template<class T, class Iter_T>
void computeStats(Iter_T first, Iter_T last, T& sum, T& median, T& mean,
                  T& var, T& std_dev, T& skew, T& kurt)
{

    size_t cnt = distance(first, last);
    sum = accumulate(first, last, T( ));
    median = computeMedian(first, last, sum);
    mean = sum / cnt;
    var = computeVariance(first, last, mean);
    std_dev = sqrt(var);
    skew = computeSkew(first, last, mean);
    kurt = computeKurtosisExcess(first, last, mean);
}

template <class T>
class pc {
public:
    T pearson_correlation(vector<T> x, vector<T> y, int n);
};

template <>
class pc <bool> {
public:
    bool pearson_correlation(vector<bool> x, vector<bool> y, int n) {
        return false;
    }
};

// taken from C++ Cookbook
template <>
class pc <float> {
public:
    float pearson_correlation(vector<float> x, vector<float> y, int n) {

        float result;
        float xmean;
        float ymean;
        float s;
        float xv;
        float yv;
        float t1;
        float t2;
        int i;

        xv = 0;
        yv = 0;
        if (n<=1) {
            //cerr << "NULL" << endl;
            result = 0;
            return result;
        }

        // Mean
        xmean = 0;
        ymean = 0;
        for (i = 0; i <= n-1; i++) {
            xmean = xmean+x.at(i);
            ymean = ymean+y.at(i);
        }
        xmean = xmean/n;
        ymean = ymean/n;

        // numerator and denominator
        s = 0;
        xv = 0;
        yv = 0;
        for (i = 0; i <= n-1; i++) {
            t1 = x.at(i)-xmean;
            t2 = y.at(i)-ymean;
            xv = xv+t1*t1;
            yv = yv+t2*t2;
            s = s+t1*t2;
        }
        if ( xv==0||yv==0 ) {
            result = 0;
        }
        else {
            result = s/(sqrt(xv)*sqrt(yv));
        }
        return result;
    }
};

#endif
