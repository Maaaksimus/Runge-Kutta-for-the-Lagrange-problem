#include "mylib.h"

int main() {
	FILE *out, *func;
	vector<double> u;
	int n = 20, max_rate = 0;
	char *endptr;
	double h, err_rate = 0, solve_rate = 0, e[2];

	printf("Enter the number of points: ");
	scanf("%d", &n);

	h = 1. / (n - 1);
	
	// shooting(n, e);
	RungeKutta();

	out = fopen("res.txt", "w");
	func = fopen("targ.txt", "w");

	max_rate = 0;

	for (int i = 0; i < n; i ++) {
		fprintf(out, "%lf %lf\n", i*h, u[0][i]);
	}

	fclose(out);
	fclose(func);

	printf("Error rate: %le\n", sqrt(err_rate));
	printf("Check results in the file\n");
	return 0;
}

void RungeKutta() {

}