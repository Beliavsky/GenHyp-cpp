#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <boost/math/special_functions/bessel.hpp>

/*
 * Program: Generalized Hyperbolic Distribution Validator
 *
 * Description:
 *   Computes the probability density at x=0, mean, and variance of the generalized
 *   hyperbolic distribution using various parameter sets. It uses the Boost.Math library's
 *   modified Bessel functions to evaluate the theoretical formulas, and also computes numerical
 *   integrals (via Simpson's rule) to verify normalization, mean, and variance.
 *
 *   The theoretical variance is computed using the variance-mean mixture representation:
 *       Var[X] = E[V] + beta^2 * Var[V],
 *   where
 *       E[V] = (delta/gamma) * K_{lambda+1}(delta*gamma)/K_lambda(delta*gamma),
 *       Var[V] = (delta^2/gamma^2) * [K_{lambda+2}(delta*gamma)/K_lambda(delta*gamma) - (K_{lambda+1}(delta*gamma)/K_lambda(delta*gamma))^2],
 *       gamma = sqrt(alpha^2 - beta^2).
 *
 * Usage:
 *   Modify the parameter sets below as needed. The program prints the theoretical and
 *   numerical values for each parameter set.
 */

// Computes the density of the generalized hyperbolic distribution at x.
double gh_density(double x, double lambda, double alpha, double beta, double delta, double mu) {
	using boost::math::cyl_bessel_k;
	const double pi = 3.14159265358979323846;
	double gamma = std::sqrt(alpha * alpha - beta * beta);
	double r = std::sqrt(delta * delta + (x - mu) * (x - mu));

    // Constant term: (gamma/delta)^lambda / (sqrt(2*pi)*K_lambda(delta*gamma))
	double A = std::pow(gamma / delta, lambda) / (std::sqrt(2 * pi) * cyl_bessel_k(lambda, delta * gamma));
    // Adjustment factor: (alpha/r)^(0.5 - lambda)
	double B = std::pow(alpha / r, 0.5 - lambda);

	double density = A * B * cyl_bessel_k(lambda - 0.5, alpha * r) * std::exp(beta * (x - mu));
	return density;
}

// Computes the theoretical mean of the GH distribution.
double gh_mean(double lambda, double alpha, double beta, double delta, double mu) {
	using boost::math::cyl_bessel_k;
	double gamma = std::sqrt(alpha * alpha - beta * beta);
	double ratio = cyl_bessel_k(lambda + 1, delta * gamma) / cyl_bessel_k(lambda, delta * gamma);
	return mu + (delta * beta / gamma) * ratio;
}

// Computes the theoretical variance of the GH distribution.
double gh_variance(double lambda, double alpha, double beta, double delta, double mu) {
	using boost::math::cyl_bessel_k;
	double gamma = std::sqrt(alpha * alpha - beta * beta);
	double ratio1 = cyl_bessel_k(lambda + 1, delta * gamma) / cyl_bessel_k(lambda, delta * gamma);
	double ratio2 = cyl_bessel_k(lambda + 2, delta * gamma) / cyl_bessel_k(lambda, delta * gamma);
	double EV = delta / gamma * ratio1;
	double VarV = (delta * delta / (gamma * gamma)) * (ratio2 - ratio1 * ratio1);
	return EV + beta * beta * VarV;
}

// Composite Simpson's rule integration for function f over [a,b] using n subintervals.
double simpson_integration(std::function<double(double)> f, double a, double b, int n) {
	if (n % 2 == 1) n++;  // Simpson's rule requires an even number of subintervals
	double h = (b - a) / n;
	double s = f(a) + f(b);
	for (int i = 1; i < n; ++i) {
		double x = a + i * h;
		s += (i % 2 == 1 ? 4 : 2) * f(x);
	}
	return s * h / 3;
}

// Structure to hold a parameter set.
struct GHParams {
	double lambda;
	double alpha;
	double beta;
	double delta;
	double mu;
};

int main() {
    // Define several parameter sets to test.
	std::vector<GHParams> paramsList = {
		{1.0, 2.0, 0.5, 1.0, 0.0},
		{0.5, 3.0, 0.8, 2.0, 1.0},
		{1.5, 2.5, 0.3, 1.5, -0.5},
		{2.0, 4.0, 1.0, 2.0, 0.0}
	};

    // Integration settings
	double L = 10.0; // integration interval: [mu - L, mu + L]
	int n_subintervals = 1000; // must be even for Simpson's rule

    // Loop over each parameter set.
	for (size_t i = 0; i < paramsList.size(); ++i) {
		const auto& p = paramsList[i];
		double a = p.mu - L;
		double b = p.mu + L;

		std::cout << "--------------------------------------------\n";
		std::cout << "Parameter Set " << i+1 << ":\n";
		std::cout << "  lambda = " << p.lambda << "\n";
		std::cout << "  alpha  = " << p.alpha  << "\n";
		std::cout << "  beta   = " << p.beta   << "\n";
		std::cout << "  delta  = " << p.delta  << "\n";
		std::cout << "  mu     = " << p.mu     << "\n\n";

		double density = gh_density(0.0, p.lambda, p.alpha, p.beta, p.delta, p.mu);
		double theo_mean = gh_mean(p.lambda, p.alpha, p.beta, p.delta, p.mu);
		double theo_variance = gh_variance(p.lambda, p.alpha, p.beta, p.delta, p.mu);

		std::cout << "Theoretical Results (at x=0 for density):\n";
		std::cout << "  Density at x = 0: " << density << "\n";
		std::cout << "  Mean     = " << theo_mean << "\n";
		std::cout << "  Variance = " << theo_variance << "\n";

	// Numerical integration to verify the normalization, mean and variance.
		auto f_norm = [&](double x) { return gh_density(x, p.lambda, p.alpha, p.beta, p.delta, p.mu); };
		double norm_integral = simpson_integration(f_norm, a, b, n_subintervals);

		auto f_mean = [&](double x) { return x * gh_density(x, p.lambda, p.alpha, p.beta, p.delta, p.mu); };
		double int_mean = simpson_integration(f_mean, a, b, n_subintervals);

		auto f_variance = [&](double x) { return (x - int_mean) * (x - int_mean) * gh_density(x, p.lambda, p.alpha, p.beta, p.delta, p.mu); };
		double int_variance = simpson_integration(f_variance, a, b, n_subintervals);

		std::cout << "\nNumerical Integration Results (over [" << a << ", " << b << "]):\n";
		std::cout << "  Integral of density (should be ~1): " << norm_integral << "\n";
		std::cout << "  Numerical Mean  = " << int_mean << "\n";
		std::cout << "  Numerical Variance = " << int_variance << "\n";
		std::cout << "--------------------------------------------\n\n";
	}

	return 0;
}
