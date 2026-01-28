#include "roots.hpp"
#include <cmath>

bool bisection(std::function<double(double)> f,
                double a, double b,
                double *root) 
{
    const double tolerance = 1e-6;
    const int max_iterations = 1000000;

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) return false; // verifies that fa and fb have opposite signs

    for (int i = 0; i < max_iterations; ++i) // calculation loop
    {
        double c = 0.5 * (a + b); // calculating c
        double fc = f(c); // calculating f(c)

        if (std::fabs(fc) <= tolerance || std::fabs(b-a) <= tolerance) // when f(c) is within the tolerance it deems c a root
        {
            *root = c;
            return true;
        }

        else if (fa * fc < 0.0) // confirms that fa and fc have opposite signs
        {
            b = c;
            fb = fc;
        }
        else 
        {
            a = c;
            fa = fc;
        }
    }

    *root = 0.5 * (a + b); // if iterations is exceeded, declare the closest values found

    return true;

}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    const double tolerance = 1e-6;
    const int max_iterations = 1000000;

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0 ) return false;
    
    for (int i = 0; i < max_iterations; ++i) 
    {
        double c = a - (fa * (b-a)) / (fb - fa); // regula falsi equation
        double fc = f(c);

        if (std::fabs(fc) <= tolerance)
        {
            *root = c;
            return true;
        }

        if (fc * fa < 0.0)
        {
            b = c;
            fb = fc;
        }

        else
        {
            a = c;
            fa = fc;
        }

    }

    *root = a - (f(a) * (b - a)) / (fb - fa);
    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    const double tolerance = 1e-6;
    const int max_iterations = 1000000;

    double x = c;

    for (int i = 0; i < max_iterations; ++i)
    {
        double fx = f(x); // calculate f(x)
        double fx_deriv = g(x); // calculate f'(x)

        if(fx_deriv == 0.0) return false; // case fails when deriv = 0

        double x_new = x - fx / fx_deriv; // calculate x_new

        if (x_new < a || x_new > b) return false; // checks if x_new stays within bounds

        if (std::fabs(x_new - x) <= tolerance) // checks if the difference between x_new and x are within the tolerance, meaning x_new is close enough
        {
            *root = x_new;
            return true;
        }

        x = x_new; // updates x to increase iteration of x

    }

    *root = x; // outputs root if max iterations exceeded
    return true;

}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    const float tolerance = 1e-6;
    const int max_iterations = 1000000;

    double x0 = c;
    double x1 = b;
    
    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < max_iterations; ++i)
    {
        if (f1 - f0 == 0) return false;

        double x2 = x1 - ((x0 - x1) / (f0 - f1)) * f1; // equation

        if (x2 < a || x2 > b) return false; // if x2 leaves bounds return that root was not found

        if(std::fabs(x2 - x1) <= tolerance) // root within tolerance, so output value
        {
            *root = x2;
            return true;
        }

        x0 = x1; // update previous values for next iteration
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    *root = x1; // output root if iterations exceeded don't get to desired tolerance
    return true;
}